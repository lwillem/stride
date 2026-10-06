#!/bin/bash
# PreToolUse hook (Bash): builds and test runs belong in the test-runner subagent.
# Denies them in the main session; tool calls from inside a subagent carry `agent_id`
# in the hook input and pass through.

input=$(cat)
[ -n "$(echo "$input" | jq -r '.agent_id // empty')" ] && exit 0

cmd=$(echo "$input" | jq -r '.tool_input.command // empty')
[ -n "$cmd" ] || exit 0

# Ignore text that is data, not commands: heredoc bodies (e.g. commit messages) and
# quoted strings. Without this, a commit message mentioning "...; gtester 22/22" is denied.
cmd=$(printf '%s\n' "$cmd" | perl -0pe '
  s/<<-?\s*(["\x27]?)(\w+)\1[^\n]*\n.*?\n\s*\2[ \t]*(\n|$)/<<HEREDOC\n/gs;
  s/\x27[^\x27]*\x27/Q/g;
  s/"(?:[^"\\]|\\.)*"/Q/g;
')

# Only match at the start of a command segment, so e.g. `grep gtester file` still works.
start='(^|[;&|(]|\bnohup|\btime|\bcd [^;&|]*&&)[[:space:]]*'
patterns=(
  "${start}([^[:space:];&|]*/)?(gtester|ctest)([[:space:];&|)]|$)"
  "${start}Rscript[[:space:]]+[^;&|]*rStride_gtester"
  "${start}make([[:space:]]+-[^[:space:]]+)*[[:space:]]+(configure|all|install|libstride|stride|gtester)([[:space:];&|)]|$)"
  "${start}make([[:space:]]+-[^[:space:]]+)*[[:space:]]*([;&|)]|$)"
  "${start}cmake[[:space:]]+--(build|install)"
  "${start}(bash[[:space:]]+)?([^[:space:];&|]*/)?stride-test\.sh"
)
for p in "${patterns[@]}"; do
  if echo "$cmd" | grep -qE "$p"; then
    jq -n '{hookSpecificOutput: {hookEventName: "PreToolUse", permissionDecision: "deny",
      permissionDecisionReason: "Builds and test runs go through the test-runner subagent (see CLAUDE.md). Delegate this to test-runner with what to run (gtester, rstride or both), an optional gtest filter and an optional install prefix."}}'
    exit 0
  fi
done
exit 0
