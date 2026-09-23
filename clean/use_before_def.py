#!/usr/bin/env python3
"""USE-BEFORE-DEFINITION check for option variables.

    usage: use_before_def.py <study script dir>

A guarded default protects an option only if it runs BEFORE the option is read.
Insert the guard beside one use and miss an earlier one, and the script dies on
"Unrecognized function or variable" - at the earlier line, which may be near the
END of a long script, after all the expensive work has been done and before
anything is saved.

That is not hypothetical: contrast_objects_tag was guarded beside
savefilenamedata while a printhdr eight lines above also read it, and prep_3
died there after 75 minutes of ComBat and contrast formation, saving nothing.

checkcode does not catch it - the line is syntactically perfect, and MATLAB
cannot know whether an undefined name will exist at run time.

Reported per option: first definition (assignment or ~exist guard) against
first bare read, in each script.
"""
import pathlib, re, sys
from collections import defaultdict

D = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else '.')

def strip_comment(line):
    out, q = [], False
    for ch in line:
        if ch == "'": q = not q
        if ch == '%' and not q: break
        out.append(ch)
    return ''.join(out)

# option-ish names: lowercase_with_underscores, assigned somewhere in the dir
NAME = re.compile(r'\b([a-z][a-z0-9_]{3,})\b')
rows = []
for f in sorted(D.glob('*.m')):
    src = [strip_comment(l) for l in f.read_text(errors='replace').splitlines()]
    assigned = set()
    for l in src:
        m = re.match(r'\s*([a-z][a-z0-9_]{3,})\s*=\s*[^=]', l)
        if m: assigned.add(m.group(1))
    guards = defaultdict(list); defs = defaultdict(list); uses = defaultdict(list)
    for i, l in enumerate(src, 1):
        for g in re.findall(r"~exist\('([a-z][a-z0-9_]{3,})'", l):
            guards[g].append(i); defs[g].append(i)
        m = re.match(r'\s*([a-z][a-z0-9_]{3,})\s*=\s*[^=]', l)
        if m: defs[m.group(1)].append(i)
        for n in NAME.findall(l):
            if n in assigned and (not m or n != m.group(1)) and "~exist('%s'" % n not in l:
                uses[n].append(i)
    for opt in sorted(set(guards)):        # only options that rely on a guard
        u = uses.get(opt, []); d = defs.get(opt, [])
        if not u or not d: continue
        if min(u) < min(d):
            rows.append((f.name, opt, min(d), min(u)))

print('  USE-BEFORE-DEFINITION check: %d script(s) in %s\n' % (len(list(D.glob('*.m'))), D))
if not rows:
    print('  PASS - every guarded option is defined before it is first read')
else:
    print('  %-52s %-26s %8s %8s' % ('SCRIPT', 'OPTION', 'DEF L', 'USE L'))
    print('  ' + '-'*98)
    for f, o, d, u in rows:
        print('  %-52s %-26s %8d %8d' % (f[:52], o, d, u))
    print('\n  %d case(s). The script will die at the USE line with' % len(rows))
    print('  "Unrecognized function or variable" - possibly after hours of work.')
    sys.exit(1)
