#!/usr/bin/env python3
"""Analyze MCP tool usage logs from Claude Code hooks.

Reads JSONL files produced by the PostToolUse logging hook and generates
summary statistics about tool usage patterns, chain sequences, failure
modes, and parameter distributions.

Usage:
    python scripts/analyze_usage.py ~/.claude/logs/multiomics-kg-usage.jsonl
    python scripts/analyze_usage.py ~/.claude/logs/multiomics-kg-usage.jsonl --output report.md
"""

import argparse
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


def load_events(path: Path) -> list[dict]:
    """Load JSONL events, skipping malformed lines."""
    events = []
    with open(path) as f:
        for i, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            try:
                events.append(json.loads(line))
            except json.JSONDecodeError:
                print(f"Warning: skipping malformed line {i}", file=sys.stderr)
    return events


def group_by_session(events: list[dict]) -> dict[str, list[dict]]:
    """Group events by session_id, preserving order."""
    sessions: dict[str, list[dict]] = defaultdict(list)
    for event in events:
        sid = event.get("session_id", "unknown")
        sessions[sid].append(event)
    return dict(sessions)


def extract_chains(session_events: list[dict]) -> list[list[str]]:
    """Extract tool call chains from a session.

    A chain is a sequence of tool calls. Chains are split when there's
    a gap of more than 5 minutes between calls (suggesting a new prompt).
    """
    chains: list[list[str]] = []
    current_chain: list[str] = []
    last_time = None

    for event in session_events:
        tool = event.get("tool_name", "unknown")
        ts = event.get("timestamp")

        if ts and last_time:
            try:
                current = datetime.fromisoformat(ts.replace("Z", "+00:00"))
                previous = datetime.fromisoformat(last_time.replace("Z", "+00:00"))
                gap = (current - previous).total_seconds()
                if gap > 300:  # 5-minute gap = new chain
                    if current_chain:
                        chains.append(current_chain)
                    current_chain = []
            except (ValueError, TypeError):
                pass

        # Strip the mcp__multiomics-kg__ prefix for readability
        short_name = tool.replace("mcp__multiomics-kg__", "")
        current_chain.append(short_name)
        last_time = ts

    if current_chain:
        chains.append(current_chain)

    return chains


def analyze(events: list[dict]) -> dict:
    """Produce analysis summary from events."""
    sessions = group_by_session(events)

    # Tool frequency
    tool_counts = Counter(
        e.get("tool_name", "unknown").replace("mcp__multiomics-kg__", "")
        for e in events
    )

    # Failure analysis
    failures = [e for e in events if e.get("event") == "PostToolUseFailure"]
    failure_counts = Counter(
        e.get("tool_name", "unknown").replace("mcp__multiomics-kg__", "")
        for e in failures
    )

    # Truncation tracking
    truncated_calls = [e for e in events if e.get("truncated") is True]
    truncation_rate = len(truncated_calls) / len(events) if events else 0

    # run_cypher usage (escape hatch signal)
    cypher_count = tool_counts.get("run_cypher", 0)
    cypher_rate = cypher_count / len(events) if events else 0

    # Chain analysis
    all_chains: list[list[str]] = []
    for session_events in sessions.values():
        all_chains.extend(extract_chains(session_events))

    # Common chain patterns (bigrams)
    bigrams = Counter()
    for chain in all_chains:
        for i in range(len(chain) - 1):
            bigrams[(chain[i], chain[i + 1])] += 1

    # Parameter distributions (what organisms/experiments are queried most)
    organisms = Counter()
    for e in events:
        inp = e.get("tool_input", {})
        if isinstance(inp, dict):
            org = inp.get("organism") or inp.get("organisms")
            if isinstance(org, str):
                organisms[org] += 1
            elif isinstance(org, list):
                for o in org:
                    organisms[o] += 1

    return {
        "total_events": len(events),
        "total_sessions": len(sessions),
        "total_failures": len(failures),
        "total_chains": len(all_chains),
        "tool_counts": dict(tool_counts.most_common()),
        "failure_counts": dict(failure_counts.most_common()),
        "truncation_rate": round(truncation_rate, 3),
        "cypher_escape_rate": round(cypher_rate, 3),
        "top_bigrams": [
            {"from": a, "to": b, "count": c}
            for (a, b), c in bigrams.most_common(15)
        ],
        "organism_distribution": dict(organisms.most_common()),
        "avg_chain_length": (
            round(sum(len(c) for c in all_chains) / len(all_chains), 1)
            if all_chains
            else 0
        ),
    }


def format_markdown(analysis: dict) -> str:
    """Format analysis as a markdown report."""
    lines = [
        "# MCP Tool Usage Report",
        "",
        f"**Total events:** {analysis['total_events']}",
        f"**Sessions:** {analysis['total_sessions']}",
        f"**Failures:** {analysis['total_failures']}",
        f"**Tool chains:** {analysis['total_chains']} "
        f"(avg length: {analysis['avg_chain_length']})",
        "",
        "## Tool Usage",
        "",
        "| Tool | Count |",
        "|------|-------|",
    ]
    for tool, count in analysis["tool_counts"].items():
        lines.append(f"| {tool} | {count} |")

    lines.extend([
        "",
        "## Key Metrics",
        "",
        f"- **Truncation rate:** {analysis['truncation_rate']:.1%} "
        f"of calls return truncated results",
        f"- **run_cypher escape rate:** {analysis['cypher_escape_rate']:.1%} "
        f"of calls use raw Cypher instead of structured tools",
    ])

    if analysis["failure_counts"]:
        lines.extend([
            "",
            "## Failures by Tool",
            "",
            "| Tool | Failures |",
            "|------|----------|",
        ])
        for tool, count in analysis["failure_counts"].items():
            lines.append(f"| {tool} | {count} |")

    if analysis["top_bigrams"]:
        lines.extend([
            "",
            "## Common Tool Transitions",
            "",
            "| From | To | Count |",
            "|------|----|-------|",
        ])
        for bg in analysis["top_bigrams"]:
            lines.append(f"| {bg['from']} | {bg['to']} | {bg['count']} |")

    if analysis["organism_distribution"]:
        lines.extend([
            "",
            "## Organism Distribution",
            "",
            "| Organism | Queries |",
            "|----------|---------|",
        ])
        for org, count in analysis["organism_distribution"].items():
            lines.append(f"| {org} | {count} |")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(
        description="Analyze MCP tool usage logs"
    )
    parser.add_argument("logfile", type=Path, help="Path to JSONL log file")
    parser.add_argument(
        "--output", "-o", type=Path, help="Output markdown report path"
    )
    parser.add_argument(
        "--json", action="store_true", help="Output raw JSON instead of markdown"
    )
    args = parser.parse_args()

    if not args.logfile.exists():
        print(f"Log file not found: {args.logfile}", file=sys.stderr)
        print("No usage data collected yet. Use Claude Code with the "
              "multiomics-research plugin to start logging.", file=sys.stderr)
        sys.exit(1)

    events = load_events(args.logfile)
    if not events:
        print("No events found in log file.", file=sys.stderr)
        sys.exit(1)

    analysis = analyze(events)

    if args.json:
        output = json.dumps(analysis, indent=2)
    else:
        output = format_markdown(analysis)

    if args.output:
        args.output.write_text(output)
        print(f"Report written to {args.output}")
    else:
        print(output)


if __name__ == "__main__":
    main()
