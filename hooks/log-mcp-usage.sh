#!/usr/bin/env bash
# Log MCP tool usage to JSONL for analysis.
# Called by Claude Code hooks on PostToolUse and PostToolUseFailure events.
# Receives hook event JSON on stdin.
#
# Output: appends one JSON line per event to ~/.claude/logs/multiomics-kg-usage.jsonl
# Fields: timestamp, event, session_id, tool_name, tool_input, response_keys,
#          truncated, total_matching, returned, error (if failure)

set -euo pipefail

LOG_DIR="${HOME}/.claude/logs"
LOG_FILE="${LOG_DIR}/multiomics-kg-usage.jsonl"
mkdir -p "$LOG_DIR"

# Read full stdin (hook event JSON)
INPUT=$(cat)

# Extract fields with jq, adding timestamp
echo "$INPUT" | jq -c '{
  timestamp: (now | todate),
  event: .hook_event_name,
  session_id: .session_id,
  tool_name: .tool_name,
  tool_input: .tool_input,
  response_keys: (.tool_response | if type == "object" then keys else null end),
  truncated: (.tool_response.truncated // null),
  total_matching: (.tool_response.total_matching // null),
  returned: (.tool_response.returned // null),
  error: (.error // null)
}' >> "$LOG_FILE"

# Exit 0 = success, non-blocking
exit 0
