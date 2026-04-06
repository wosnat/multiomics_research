#!/usr/bin/env python3
"""Run trigger evaluation for the research skill.

Tests whether the research skill triggers (is invoked by Claude) for
queries where it should, and does NOT trigger for queries where it shouldn't.

Adapted from the Anthropic skill-creator pattern:
  https://github.com/anthropics/skills/tree/main/skills/skill-creator

Usage:
    python scripts/run_trigger_eval.py evals/trigger_evals.json --skill-path skills/research-methodology
    python scripts/run_trigger_eval.py evals/trigger_evals.json --skill-path skills/research-methodology --runs-per-query 3
"""

import argparse
import json
import os
import subprocess
import sys
import uuid
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


def find_project_root() -> Path:
    """Find the project root by walking up from cwd looking for .claude/."""
    current = Path.cwd()
    for parent in [current, *current.parents]:
        if (parent / ".claude").is_dir():
            return parent
    return current


def run_single_query(
    query: str,
    skill_name: str,
    skill_description: str,
    timeout: int,
    project_root: str,
) -> bool:
    """Run a single query and return whether the skill was triggered.

    Creates a temporary command file, runs `claude -p`, and checks
    if Claude invoked the Skill tool for our skill.
    """
    unique_id = uuid.uuid4().hex[:8]
    clean_name = f"{skill_name}-trigger-{unique_id}"
    project_commands_dir = Path(project_root) / ".claude" / "commands"
    command_file = project_commands_dir / f"{clean_name}.md"

    try:
        project_commands_dir.mkdir(parents=True, exist_ok=True)
        indented_desc = "\n  ".join(skill_description.split("\n"))
        command_content = (
            f"---\n"
            f"description: |\n"
            f"  {indented_desc}\n"
            f"---\n\n"
            f"# {skill_name}\n\n"
            f"This skill handles: {skill_description}\n"
        )
        command_file.write_text(command_content)

        cmd = [
            "claude", "-p", query,
            "--output-format", "stream-json",
            "--verbose",
        ]

        env = {k: v for k, v in os.environ.items() if k != "CLAUDECODE"}

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=project_root,
            env=env,
        )

        # Check if the skill was triggered in the output
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
                    if content.get("type") != "tool_use":
                        continue
                    tool_name = content.get("name", "")
                    tool_input = content.get("input", {})
                    if tool_name == "Skill" and clean_name in tool_input.get("skill", ""):
                        return True
                    if tool_name == "Read" and clean_name in tool_input.get("file_path", ""):
                        return True

        return False

    except subprocess.TimeoutExpired:
        return False
    finally:
        if command_file.exists():
            command_file.unlink()


def parse_skill_md(skill_path: Path) -> tuple[str, str]:
    """Extract name and description from SKILL.md frontmatter."""
    skill_file = skill_path / "SKILL.md"
    text = skill_file.read_text()

    name = ""
    description = ""
    in_frontmatter = False

    for line in text.splitlines():
        if line.strip() == "---":
            if in_frontmatter:
                break
            in_frontmatter = True
            continue
        if in_frontmatter:
            if line.startswith("name:"):
                name = line.split(":", 1)[1].strip()
            elif line.startswith("description:"):
                description = line.split(":", 1)[1].strip()

    return name, description


def main():
    parser = argparse.ArgumentParser(description="Run trigger evaluation")
    parser.add_argument("eval_set", type=Path, help="Path to trigger eval JSON")
    parser.add_argument("--skill-path", type=Path, required=True, help="Path to skill directory")
    parser.add_argument("--num-workers", type=int, default=5, help="Parallel workers")
    parser.add_argument("--timeout", type=int, default=30, help="Timeout per query (seconds)")
    parser.add_argument("--runs-per-query", type=int, default=1, help="Runs per query")
    parser.add_argument("--trigger-threshold", type=float, default=0.5, help="Trigger rate threshold")
    args = parser.parse_args()

    eval_set = json.loads(args.eval_set.read_text())
    name, description = parse_skill_md(args.skill_path)
    project_root = find_project_root()

    print(f"Evaluating skill: {name}", file=sys.stderr)
    print(f"Description: {description}", file=sys.stderr)
    print(f"Queries: {len(eval_set)} ({args.runs_per_query} runs each)", file=sys.stderr)
    print("", file=sys.stderr)

    results = []
    with ProcessPoolExecutor(max_workers=args.num_workers) as executor:
        future_to_info = {}
        for item in eval_set:
            for run_idx in range(args.runs_per_query):
                future = executor.submit(
                    run_single_query,
                    item["query"],
                    name,
                    description,
                    args.timeout,
                    str(project_root),
                )
                future_to_info[future] = (item, run_idx)

        query_triggers: dict[str, list[bool]] = {}
        query_items: dict[str, dict] = {}
        for future in as_completed(future_to_info):
            item, _ = future_to_info[future]
            query = item["query"]
            query_items[query] = item
            if query not in query_triggers:
                query_triggers[query] = []
            try:
                query_triggers[query].append(future.result())
            except Exception as e:
                print(f"Warning: query failed: {e}", file=sys.stderr)
                query_triggers[query].append(False)

    for query, triggers in query_triggers.items():
        item = query_items[query]
        trigger_rate = sum(triggers) / len(triggers)
        should_trigger = item["should_trigger"]
        if should_trigger:
            did_pass = trigger_rate >= args.trigger_threshold
        else:
            did_pass = trigger_rate < args.trigger_threshold
        results.append({
            "query": query,
            "should_trigger": should_trigger,
            "trigger_rate": trigger_rate,
            "triggers": sum(triggers),
            "runs": len(triggers),
            "pass": did_pass,
            "reason": item.get("reason", ""),
        })

    passed = sum(1 for r in results if r["pass"])
    total = len(results)

    print(f"Results: {passed}/{total} passed", file=sys.stderr)
    for r in results:
        status = "PASS" if r["pass"] else "FAIL"
        rate = f"{r['triggers']}/{r['runs']}"
        expected = "should" if r["should_trigger"] else "should NOT"
        print(f"  [{status}] rate={rate} ({expected} trigger): {r['query'][:60]}", file=sys.stderr)

    output = {
        "skill_name": name,
        "description": description,
        "summary": {"total": total, "passed": passed, "failed": total - passed},
        "results": results,
    }
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
