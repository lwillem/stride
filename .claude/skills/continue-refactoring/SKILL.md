---
name: continue-refactoring
description: Resume the rStride/STRIDE refactoring from the plan's "State of play" section, propose the next step, and keep the plan updated. Use when the user says continue/resume the refactoring, or runs /continue-refactoring.
---

# Continue the refactoring

The plan is `doc/markdown/rStride_refactoring_plan.md` (~2000 lines). Read it selectively.

## 1. Get oriented (read only what's needed)
1. Read **section 0 "State of play"** only: from the `## 0.` heading up to `## 1.`
   (find the line numbers with `grep -n "^## " doc/markdown/rStride_refactoring_plan.md`).
2. Check that the repo matches it: `git branch --show-current`, `git status --short`,
   `git log --oneline -5`. Point out any mismatch (other branch, uncommitted work, commits
   that section 0 doesn't mention).
3. For the first item under "Next, in order", read **only** its phase section in section 2
   and the findings (F-numbers) it cites in section 1. Look up structures in
   `rStride_architecture.md` by heading only when needed.

## 2. Propose before acting
Give the user a short briefing (≤15 lines):
- where things stand (branch, last completed step)
- the proposed next step, its scope, and the commits you expect to make
- any open question in section 0 that blocks it

Wait for the user to confirm or pick another item before you change any code.

## 3. While working
- Follow the plan's **Working rules** (end of section 2): never mix a behaviour change and a
  structure change in one commit; every phase ends green.
- Run every build and test through the `test-runner` subagent, never in this session.
- Commit only when the user asks.
- Never commit or push to `master`. Work on a branch (create one before the first commit),
  push it and open a pull request with `gh pr create`; merge only when the user asks.

## 4. Keep the plan current with every commit
Every refactoring commit updates the plan in the same commit, or in a `docs:` commit
right after it, so the plan is never more than one step behind however the session ends:
- **Section 0 "State of play"**: the date, branches, "Done, with evidence", "Next, in order"
  and open questions.
- **Phase status** in its section-2 heading (e.g. `— COMPLETE`, `— steps 1-2 done`).

A Stop hook (`.claude/hooks/check-state-of-play.sh`) checks this: if HEAD has commits newer
than the last change to the plan, it asks once for the update. If those commits are
unrelated to the refactoring, say so in one line instead.

## 5. Architecture document (optional, ask first)
If a step changed a structure that `rStride_architecture.md` describes (install flow,
experiment pipeline, contact pools, transmission/calibration, build, test suites, repo
layout), name the affected section and ask the user whether to update it now, later, or not
at all. Don't edit it without a yes.

## 6. Before stopping (or when the context gets large)
Check that section 0 matches the last commit, then suggest starting a new session rather
than continuing in a long one.
