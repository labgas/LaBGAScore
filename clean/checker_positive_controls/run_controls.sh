#!/bin/bash
# Verify both script checkers still catch what they were written to catch.
#
# Each control is POSITIVE for its own checker and NEGATIVE for the other, so a
# correct run reports the matching file(s) per checker and nothing else.
# Both checkers exit non-zero when they find something, which here is SUCCESS.
#
# use_before_def.py has TWO controls, one per guard spelling: the three-line
# "if ~exist(...)\n x = ...\nend" and the one-line
# "if ~exist('x','var'), x = []; end". Only the first was covered until
# 2026-10-05, and the second is the form the real templates actually use - the
# checker passed all eleven guarded options in prep_3a without examining any of
# them. Keep both.
#
#   usage: clean/checker_positive_controls/run_controls.sh
CD="$(cd "$(dirname "$0")" && pwd)"
CLEAN="$(dirname "$CD")"
fail=0

echo "=== use_before_def.py ==="
out=$(python3 "$CLEAN/use_before_def.py" "$CD" 2>&1); echo "$out"
grep -q "control_use_before_def.m *contrast_objects_tag" <<<"$out" \
  || { echo "!! FAIL: did not flag control_use_before_def.m"; fail=1; }
grep -q "control_use_before_def_guard.m *cv_seed_mvpa_reg_cov" <<<"$out" \
  || { echo "!! FAIL: did not flag control_use_before_def_guard.m (one-line guard form)"; fail=1; }
grep -q "control_set_after_use.m" <<<"$out" \
  && { echo "!! FAIL: flagged the OTHER checker's control"; fail=1; }

echo
echo "=== set_after_use.py ==="
out=$(python3 "$CLEAN/set_after_use.py" "$CD" 2>&1); echo "$out"
grep -q "control_set_after_use.m *results_tag" <<<"$out" \
  || { echo "!! FAIL: did not flag control_set_after_use.m"; fail=1; }
grep -q "control_use_before_def" <<<"$out" \
  && { echo "!! FAIL: flagged the OTHER checker's control"; fail=1; }

echo
if [ $fail -eq 0 ]; then
  echo "BOTH CHECKERS OK - each caught its own control and ignored the other's."
else
  echo "CHECKER REGRESSION - see the !! lines above."
fi
exit $fail
