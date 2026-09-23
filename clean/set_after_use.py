#!/usr/bin/env python3
"""SET-AFTER-USE check: an option SET below the line that already consumed it.

    usage: set_after_use.py <script dir>

Distinct from use-before-definition and invisible to it. The variable IS
defined before use - typically by a guarded default near the top - so nothing
errors. The problem is that the author's real setting sits BELOW the consumer,
so the value in force is the default, and the setting silently does nothing.

The case that motivated this: the decoding template builds tdt_resultsdir from
results_tag at the '% OUTPUT DIRECTORIES' block, a few lines after

    if ~exist('results_tag','var'), results_tag = ''; end

A study copy that set results_tag further down - beside scaled_contrast_dir -
was read too late. Eleven runs wrote to the same untagged folder and overwrote
each other's maps. Statistics were unaffected (results are always recomputed);
the saved artifacts were not.

The pattern flagged is therefore:

    definition (often a guard)  ->  USE  ->  top-level assignment

Only column-0 assignments count as "the setting": indented ones are loop
counters and accumulators, and including them buried the signal under hundreds
of false positives.

ADVISORY, not a gate. One benign pattern still trips it: a variable legitimately
reused for successive outputs - savefilename for a second file, figtitle for the
next figure, varnames re-derived after covariates are excluded. Expect a handful
of those per model and read them rather than chasing them. It is worth keeping
because the failure it does catch is silent, expensive and otherwise invisible.
"""
import pathlib, re, sys
from collections import defaultdict

D = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else '.')

def strip(line):
    out, q = [], False
    for ch in line:
        if ch == "'": q = not q
        if ch == '%' and not q: break
        out.append(ch)
    return ''.join(out)

ANY_ASSIGN = re.compile(r'\s*([a-z][a-z0-9_]{3,})\s*=\s*[^=]')
TOP_ASSIGN = re.compile(r'([a-z][a-z0-9_]{3,})\s*=\s*[^=]')
GUARD      = re.compile(r"~exist\('([a-z][a-z0-9_]{3,})'")

rows = []
for f in sorted(D.glob('*.m')):
    src = [strip(l) for l in f.read_text(errors='replace').splitlines()]
    defs, uses, sets_ = defaultdict(list), defaultdict(list), defaultdict(list)
    for i, l in enumerate(src, 1):
        for g in GUARD.findall(l):
            defs[g].append(i)
        m = ANY_ASSIGN.match(l)
        if m:
            defs[m.group(1)].append(i)
            if TOP_ASSIGN.match(l):
                sets_[m.group(1)].append(i)
        for n in re.findall(r'\b([a-z][a-z0-9_]{3,})\b', l):
            if not m or n != m.group(1):
                uses[n].append(i)
    for opt, sw in sets_.items():
        d, u = defs.get(opt, []), uses.get(opt, [])
        if not d or not u:
            continue
        # a use that already happened before the last top-level setting,
        # and that use was itself covered by an earlier definition
        early_use = [x for x in u if x < max(sw) and any(y < x for y in d)]
        if early_use:
            rows.append((f.name, opt, min(d), early_use[0], max(sw)))

print('  SET-AFTER-USE check: %d script(s) in %s\n' % (len(list(D.glob('*.m'))), D))
if not rows:
    print('  PASS - no option is set after something already read it')
else:
    print('  %-44s %-20s %7s %7s %7s' % ('SCRIPT', 'OPTION', 'def L', 'USED L', 'SET L'))
    print('  ' + '-'*90)
    for f, o, d, u, s in rows:
        print('  %-44s %-20s %7d %7d %7d' % (f[:44], o, d, u, s))
    print('\n  %d case(s). The value in force at the USE line is the DEFAULT;' % len(rows))
    print('  the setting below it is silently ignored by whatever already ran.')
    sys.exit(1)
