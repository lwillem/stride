# STRIDE – notes for Claude Code

## Refactoring in progress
A multi-session refactoring is underway. The documents are in `doc/markdown/`:
- `rStride_refactoring_plan.md`: the plan. Section 0 "State of play" is the current status and the next steps.
- `rStride_architecture.md`: how the system works today (look things up there, don't read it end to end).
- Related: `immunity_clustering_plan.md`, `measles_usa_rm_discussion.md`, `measles_usa_rm_merge_result.md`, `removed_features/`.

At the start of a session, unless the first request is clearly about something else, ask:
"Continue the refactoring? I'll read the State of play in `rStride_refactoring_plan.md` first."
If yes, follow the `/continue-refactoring` skill. Don't read the full plan (~2000 lines) up front; read section 0 and then only the sections it points to.

## Builds and tests
- Always run builds and tests through the `test-runner` subagent (`.claude/agents/test-runner.md`); never run `gtester` or `rStride_gtester_covid19.R` in the main session.
- A PreToolUse hook (`.claude/hooks/guard-test-runs.sh`) denies build/test commands outside a subagent.
- Pass it what to run (`gtester`, `rstride` or `both`), an optional gtest filter, and an optional install prefix (default `~/opt/stride`). Use its summary, and open the log files it names only when a failure needs a closer look.
