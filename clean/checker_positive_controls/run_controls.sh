#!/bin/bash
# Verify both script checkers still catch what they were written to catch.
#
# Each control is POSITIVE for its own checker and NEGATIVE for the other, so a
# correct run reports exactly one case per checker, naming the matching file.
# Both checkers exit non-zero when they find something, which here is SUCCESS.
#
#   usage: clean/checker_positive_controls/run_controls.sh
CD="$(cd "$(dirname "$0")" && pwd)"
CLEAN="$(dirname "$CD")"
fail=0

echo "=== use_before_def.py ==="
out=$(python3 "$CLEAN/use_before_def.py" "$CD" 2>&1); echo "$out"
grep -q "control_use_before_def.m *contrast_objects_tag" <<<"$out" \
  || { echo "!! FAIL: did not flag control_use_before_def.m"; fail=1; }
grep -q "control_set_after_use.m" <<<"$out" \
  && { echo "!! FAIL: flagged the OTHER checker's control"; fail=1; }

echo
echo "=== set_after_use.py ==="
out=$(python3 "$CLEAN/set_after_use.py" "$CD" 2>&1); echo "$out"
grep -q "control_set_after_use.m *results_tag" <<<"$out" \
  || { echo "!! FAIL: did not flag control_set_after_use.m"; fail=1; }
grep -q "control_use_before_def.m" <<<"$out" \
  && { echo "!! FAIL: flagged the OTHER checker's control"; fail=1; }

echo
if [ $fail -eq 0 ]; then
  echo "BOTH CHECKERS OK - each caught its own control and ignored the other's."
else
  echo "CHECKER REGRESSION - see the !! lines above."
fi
exit $fail
