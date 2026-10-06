#!/bin/bash
# Stop hook: block once when the branch has commits newer than the last change to the
# refactoring plan, so section 0 "State of play" gets updated before the session stops.
# Fires at most once per HEAD (acknowledged HEADs are recorded under .git/), so declining
# an update for an unrelated commit doesn't trigger it again on every turn.

input=$(cat)
[ "$(echo "$input" | jq -r '.stop_hook_active // false')" = "true" ] && exit 0

cd "${CLAUDE_PROJECT_DIR:-.}" || exit 0
git rev-parse --git-dir >/dev/null 2>&1 || exit 0

PLAN=doc/markdown/rStride_refactoring_plan.md
[ -f "$PLAN" ] || exit 0

# Uncommitted edits to the plan count as "being updated".
git diff --quiet HEAD -- "$PLAN" 2>/dev/null || exit 0

plan_commit=$(git log -1 --format=%H -- "$PLAN")
[ -n "$plan_commit" ] || exit 0
behind=$(git rev-list --count "$plan_commit"..HEAD)
[ "$behind" -gt 0 ] || exit 0

head=$(git rev-parse HEAD)
ack="$(git rev-parse --git-dir)/claude-state-of-play-ack"
grep -qx "$head" "$ack" 2>/dev/null && exit 0
echo "$head" >> "$ack"

subjects=$(git log --format='- %h %s' "$plan_commit"..HEAD | head -10)
reason="The branch has $behind commit(s) since $PLAN was last changed:
$subjects
If any of these are refactoring work, update section 0 \"State of play\" (and the phase status in its section-2 heading) now. If they are unrelated to the refactoring, say so in one line and stop."

jq -n --arg r "$reason" '{decision: "block", reason: $r}'
