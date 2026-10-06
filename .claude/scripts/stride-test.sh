#!/bin/bash
# Build and test entry point for the test-runner subagent. One plain command per step, so a
# single permission rule covers it (compound `cd ... && ./bin/gtester > log` is not matched).
#
#   .claude/scripts/stride-test.sh build   [prefix]
#   .claude/scripts/stride-test.sh gtester [prefix] [gtest args...]
#   .claude/scripts/stride-test.sh rstride [prefix]     # starts in the background
#
# prefix defaults to $HOME/opt/stride; pass "-" to keep the default when adding gtest args.
# Full logs go to $TMPDIR; only short summaries are printed.

step=$1; shift
prefix=${1:--}; shift
[ "$prefix" = "-" ] && prefix="$HOME/opt/stride"
logdir=${TMPDIR:-/tmp}
repo=$(cd "$(dirname "$0")/../.." && pwd)

# Claude Code may itself run under Rosetta (x86_64). Its children then configure for
# x86_64, cannot link the arm64 libomp, and OpenMP silently drops out (plan F13.2,
# half the gtester suite vanishes). On Apple silicon, always run natively.
native=()
[ "$(sysctl -n hw.optional.arm64 2>/dev/null)" = "1" ] && native=(arch -arm64)

case "$step" in
  build)
    L="$logdir/stride_build.log"
    cd "$repo" || exit 1
    pfx=""
    [ "$prefix" != "$HOME/opt/stride" ] && pfx="CMAKE_INSTALL_PREFIX=$prefix"
    ("${native[@]}" make configure $pfx && "${native[@]}" make all $pfx \
      && "${native[@]}" make install $pfx) > "$L" 2>&1
    rc=$?
    echo "build exit=$rc  commit=$(git rev-parse --short HEAD)  branch=$(git branch --show-current)"
    grep -nE "error:|Error [0-9]|undefined reference|ld: " "$L" | head -15
    echo "warnings: $(grep -cE "warning:" "$L")"
    echo "log: $L"
    exit $rc
    ;;
  gtester)
    G="$logdir/stride_gtester.log"
    cd "$prefix" || exit 1
    "${native[@]}" ./bin/gtester "$@" --gtest_output=xml:tests/gtester_all.xml > "$G" 2>&1
    rc=$?
    echo "gtester exit=$rc"
    tail -5 "$G"
    grep -E "^\[  FAILED  \]" "$G" | sort -u
    echo "log: $G"
    exit $rc
    ;;
  rstride)
    R="$logdir/stride_rgtester.log"
    cd "$prefix" || exit 1
    nohup "${native[@]}" Rscript bin/rStride_gtester_covid19.R > "$R" 2>&1 &
    echo "rstride started pid=$!"
    echo "log: $R"
    ;;
  *)
    echo "usage: $0 build|gtester|rstride [prefix] [gtest args...]" >&2
    exit 2
    ;;
esac
