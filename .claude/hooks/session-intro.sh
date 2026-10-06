#!/bin/bash
# SessionStart hook (startup): show a short intro on how Claude works in this repo,
# with the current branch and the refactoring plan's next step.

cd "${CLAUDE_PROJECT_DIR:-.}" || exit 0
PLAN=doc/markdown/rStride_refactoring_plan.md

branch=$(git branch --show-current 2>/dev/null)
last=$(git log -1 --format='%h %s' 2>/dev/null)
dirty=$(git status --short 2>/dev/null | wc -l | tr -d ' ')
state=""; next=""
if [ -f "$PLAN" ]; then
  state=$(grep -m1 '^## 0\.' "$PLAN" | sed 's/^## 0\. *//')
  next=$(awk '/^### Next, in order/{f=1;next} f&&/^1\./{print;exit}' "$PLAN" | sed 's/^1\. *//; s/\*\*//g')
fi

msg="STRIDE / rStride - Claude Code setup
Intent: a multi-session refactoring of STRIDE/rStride, tracked in doc/markdown/rStride_refactoring_plan.md.
  Plan:    ${state:-not found}
  Next:    ${next:-see section 0}
  Repo:    ${branch:-?} @ ${last:-?} (${dirty} uncommitted change(s))
Skills and agents:
  /continue-refactoring  resume from section 0 'State of play' and propose the next step
  test-runner (Haiku)    runs every build, gtester and rStride regression test; the main session is blocked from running them
Hooks: on stop, a check asks for a State of play update when commits are newer than the plan.
       rStride_architecture.md is only updated after you say yes."

jq -n --arg m "$msg" '{systemMessage: $m}'
