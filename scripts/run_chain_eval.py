#!/usr/bin/env python3
"""Run tool chain evaluations against the MCP server via claude -p.

Loads chain eval cases from YAML, runs each question through Claude with
the MCP server available, captures the tool call sequence from the stream
output, and compares against expected chains and anti-patterns.

Usage:
    python scripts/run_chain_eval.py evals/chain_evals.yaml
    python scripts/run_chain_eval.py evals/chain_evals.yaml --output benchmarks/chains/
"""

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml


def run_question(question: str, timeout: int = 120) -> dict:
    """Run a question through claude -p and capture tool calls.

    Returns dict with tool_calls list and raw output.
    """
    cmd = [
        "claude", "-p", question,
        "--output-format", "stream-json",
        "--verbose",
    ]

    env = {k: v for k, v in os.environ.items() if k != "CLAUDECODE"}

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            env=env,
        )
    except subprocess.TimeoutExpired:
        return {"tool_calls": [], "error": "timeout", "raw": ""}

    # Parse stream-json output to extract tool calls
    tool_calls = []
    final_text = ""

    for line in result.stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            event = json.loads(line)
        except json.JSONDecodeError:
            continue

        if event.get("type") == "assistant":
            message = event.get("message", {})
            for content in message.get("content", []):
                if content.get("type") == "tool_use":
                    tool_calls.append({
                        "tool": content["name"],
                        "input": content.get("input", {}),
                    })
                elif content.get("type") == "text":
                    final_text += content.get("text", "")

        elif event.get("type") == "result":
            if not final_text:
                final_text = event.get("result", "")

    return {
        "tool_calls": tool_calls,
        "final_text": final_text,
        "raw_stdout_lines": len(result.stdout.splitlines()),
    }


def strip_mcp_prefix(name: str) -> str:
    """Remove mcp__multiomics-kg__ prefix."""
    return name.replace("mcp__multiomics-kg__", "")


def evaluate_chain(case: dict, result: dict) -> dict:
    """Evaluate a single chain result against expected patterns."""
    actual_tools = [
        strip_mcp_prefix(tc["tool"])
        for tc in result["tool_calls"]
        if "mcp__multiomics-kg__" in tc["tool"]
    ]
    actual_sequence = " -> ".join(actual_tools) if actual_tools else "(no MCP tools called)"

    # Check expected chain (order matters, but gaps are OK)
    expected = case.get("expected_chain", [])
    chain_matches = []
    search_from = 0
    for step in expected:
        expected_tool = step["tool"]
        found = False
        for i in range(search_from, len(actual_tools)):
            if actual_tools[i] == expected_tool:
                # Check params_contain if specified
                if "params_contain" in step:
                    actual_input = result["tool_calls"][i].get("input", {})
                    params_ok = all(
                        str(actual_input.get(k)) == str(v)
                        or (isinstance(v, list) and actual_input.get(k) == v)
                        for k, v in step["params_contain"].items()
                    )
                    if not params_ok:
                        chain_matches.append({
                            "tool": expected_tool,
                            "found": True,
                            "params_match": False,
                            "reason": step.get("reason", ""),
                        })
                        continue

                chain_matches.append({
                    "tool": expected_tool,
                    "found": True,
                    "params_match": True,
                    "reason": step.get("reason", ""),
                })
                search_from = i + 1
                found = True
                break
        if not found:
            chain_matches.append({
                "tool": expected_tool,
                "found": False,
                "params_match": False,
                "reason": step.get("reason", ""),
            })

    chain_pass = all(m["found"] for m in chain_matches)

    # Check anti-patterns
    anti_pattern_violations = []
    import re
    for ap in case.get("anti_patterns", []):
        pattern = ap["pattern"]
        full_sequence = " ".join(actual_tools)
        if re.search(pattern, full_sequence):
            anti_pattern_violations.append({
                "pattern": pattern,
                "reason": ap["reason"],
            })

    anti_pass = len(anti_pattern_violations) == 0

    return {
        "id": case["id"],
        "question": case["question"].strip(),
        "category": case.get("category", "unknown"),
        "actual_sequence": actual_sequence,
        "actual_tools": actual_tools,
        "chain_matches": chain_matches,
        "chain_pass": chain_pass,
        "anti_pattern_violations": anti_pattern_violations,
        "anti_pass": anti_pass,
        "overall_pass": chain_pass and anti_pass,
    }


def main():
    parser = argparse.ArgumentParser(description="Run tool chain evaluations")
    parser.add_argument("eval_file", type=Path, help="Path to chain eval YAML")
    parser.add_argument("--output", "-o", type=Path, help="Output directory for results")
    parser.add_argument("--timeout", type=int, default=120, help="Timeout per question (seconds)")
    parser.add_argument("--case", type=str, help="Run a single case by ID")
    args = parser.parse_args()

    cases = yaml.safe_load(args.eval_file.read_text())

    if args.case:
        cases = [c for c in cases if c["id"] == args.case]
        if not cases:
            print(f"Case '{args.case}' not found", file=sys.stderr)
            sys.exit(1)

    results = []
    for i, case in enumerate(cases):
        case_id = case["id"]
        print(f"[{i+1}/{len(cases)}] Running: {case_id}...", file=sys.stderr)

        run_result = run_question(case["question"], timeout=args.timeout)
        eval_result = evaluate_chain(case, run_result)
        results.append(eval_result)

        status = "PASS" if eval_result["overall_pass"] else "FAIL"
        print(f"  [{status}] {eval_result['actual_sequence']}", file=sys.stderr)

    # Summary
    passed = sum(1 for r in results if r["overall_pass"])
    total = len(results)
    print(f"\nResults: {passed}/{total} passed", file=sys.stderr)

    output = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "total": total,
        "passed": passed,
        "failed": total - passed,
        "results": results,
    }

    if args.output:
        args.output.mkdir(parents=True, exist_ok=True)
        out_file = args.output / f"chain_eval_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        out_file.write_text(json.dumps(output, indent=2))
        print(f"Results written to {out_file}", file=sys.stderr)
    else:
        print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
