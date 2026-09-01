# Provenance and dependency tooling

Records which version of every dependency produced a given analysis result, and
documents what each script depends on.

- [The problem](#the-problem)
- [Do this first](#do-this-first)
- [Recording provenance for new runs](#recording-provenance-for-new-runs)
- [Reconstructing provenance for past runs](#reconstructing-provenance-for-past-runs)
- [Dependency documentation](#dependency-documentation)
- [What the fields mean](#what-the-fields-mean)
- [What this cannot tell you](#what-this-cannot-tell-you)
- [Design notes](#design-notes)

## The problem

Analysis scripts are copied into a study's `code` subdataset and adapted there, so each
project owns a frozen copy. Their **dependencies are not frozen**: CanlabCore,
CanlabPrivate, the LaBGAS fork of CANlab_help_examples, canlab_single_trials and ~25
other clones under `/data/master_github_repos` are shared and keep moving.

Vendoring them per project was considered and rejected: Neuroimaging_Pattern_Masks
alone is 18 GB, CanlabCore 4.5 GB, and pinning would block the bug fixes that are
actively being made. Recording what was used costs nothing and blocks nothing.

## Do this first

```bash
bash clean/labgascore_prov_protect_reflogs.sh /data/master_github_repos
```

Git prunes reflog entries after 90 days by default, silently, during garbage
collection. **The reflog is the only record of which commit a clone actually had
checked out at a past moment** — it cannot be reconstructed afterwards, because a
commit's own author/committer date says nothing about when your machine moved to it.
Once pruned, the provenance of past analyses is gone for good.

The script sets `gc.reflogExpire=never` and `gc.reflogExpireUnreachable=never` on every
clone. It is idempotent and touches nothing else.

Then archive the reflogs periodically, so the record survives a re-clone:

```matlab
report = LaBGAScore_prov_protect_reflogs;
```

## Recording provenance for new runs

Replace the `publish` call documented in each script header:

```matlab
% instead of
publish('cfs_secondlevel_m14_s7_c2a_second_level_regression', 'outputDir', htmlsavedir)

% use
LaBGAScore_prov_publish('cfs_secondlevel_m14_s7_c2a_second_level_regression', htmlsavedir)
```

The HTML report gains a Provenance section; a MATLAB table and a `.tsv` go to the
model's `results/notes/`. No script templates need editing — deliberately, since the
templates live partly in the sibling repo and every study has its own renamed copies
that would never pick up such an edit.

### Figure size and screen size

> **Setting up your session:** see *"Set up your X2go display for publishing figures"* in
> [`LaBGAS_fMRI_analysis_workflow.md`](../LaBGAS_fMRI_analysis_workflow.md) for recommended
> X2go settings per screen size, and run `LaBGAScore_check_display` to check your own.


`LaBGAScore_prov_publish` deliberately does **not** set `maxWidth`/`maxHeight`, so figures
are written at exactly the resolution the bare `publish` call in each script header
produces. Those options look like display settings but are not: `publish` honours them by
resizing the `.png` on disk, so the detail is permanently gone. A 1600x900 figure published
with `maxWidth 1600, maxHeight 800` lands as **1422x800**; second-level brain montages are
1707x913 natively, which is exactly where that detail matters.

Fitting a small laptop is a display problem, and is solved in the report instead. Published
HTML from MATLAB has no viewport tag, no width attributes on `<img>`, and no `max-width` in
its stylesheet, so a 1707px montage scrolls a 1366px screen sideways. The wrapper injects a
short stylesheet after publish's own, which:

- scales images down to the window on narrow screens (`img { max-width:100%; height:auto }`)
  while leaving them at native size on a display wide enough to show them;
- lets long console output and wide stats tables scroll inside their own box instead of
  stretching the page.

Pass `'responsive', false` for byte-identical `publish` output, or set `maxWidth` explicitly
if you deliberately want smaller files (for e-mail, say) and accept the resolution loss.

Leaving the caps unset also matters for `plugin_set_figure_size` in CANlab_help_examples,
which sizes figures in inches so that the font-to-canvas ratio does not depend on the
client display. A `maxWidth` would resize the resulting `.png` afterwards and quietly undo
that work.

Because figure capture depends on the display, the snapshot records the session's screen
geometry and DPI **and the dimensions of every figure the report contains** (in
`PROV.figures`, appended to the `.mat` after publishing). Both appear in the report header:

    multifig · 2026-08-31T14:49:10+02:00 · gbw-s-labgas01 · MATLAB 2021a
      · screen 1718x1360 @ 133 DPI · figures up to 1684x1052 (display-limited)

`publish` captures what is on screen, so a figure requested larger than the display comes
back at display size. Any figure as wide as the screen is flagged `display-limited`, meaning
the machine rather than the script decided its size and the report will look different
elsewhere.

This is also why `plugin_set_figure_size` in CANlab_help_examples was changed to fit its
request to the display, preserving aspect ratio. Previously its fixed 16x10 inch default
needed 2128x1330 px, did not fit the LaBGAS server's 1718x1360 screen, and was captured at
1718x1254 - **aspect 1.37 instead of 1.60**, i.e. exactly the `WindowState maximized`
behaviour it exists to avoid, and silently, since `get(fh,'Position')` still reported
16x10. It now scales both dimensions by one factor, so the aspect is always honoured, and
repositions the window fully on screen (a window hanging off the edge is clamped at capture
however it was sized). On this server 16x10 in becomes 12.7x7.91 in, captured at 1684x1052,
aspect 1.601. It says so once per session rather than once per figure.

To record provenance outside `publish`, call the snapshot directly:

```matlab
PROV = LaBGAScore_prov_snapshot('scriptname', 'my_script', 'savedir', notesdir);
```

## Reconstructing provenance for past runs

```matlab
IDX = LaBGAScore_dep_build_index;
T = LaBGAScore_prov_resolve_retrospective('/data/proj_cfs', ...
        'index', IDX, ...
        'followrepos', {'CanlabCore','CanlabPrivate','CANlab_help_examples', ...
                        'canlab_single_trials','LaBGAScore'});
```

This works because two records survive independently:

1. every artifact carries its own date:
   - a report published by MATLAB has the run date (`DC.date`), the script name
     (`DC.source`), the MATLAB version, **and its own complete source** between
     `##### SOURCE BEGIN #####` markers;
   - a result `.mat` has a `Created on:` stamp in its 116-byte header, with a time of
     day, which survives git-annex and copying just as well
2. every clone's `.git/logs/HEAD` records each commit HEAD moved to, and when

**Result files matter as much as reports.** Only a few scripts per model are ever
published to HTML - the `prep_` scripts write `.mat` files instead - so resolving reports
alone leaves most of the pipeline undocumented. In `proj_cfs` that is 29 reports against
98 result files. Both are covered; pass `'artifacts', {'report'}` to restrict to reports.

Each result file is attributed to the script that wrote it in three steps:

1. **The CANlab naming convention** — `data_objects.mat` comes from `prep_2`, and so on.
2. **A search of the model's own scripts** for the filename, accepted only when exactly
   one mentions it. This picks up study-specific outputs — The Decoding Toolbox's
   `res_*.mat`, JuSpace's `neurotransmitter_*.mat` — with no convention table at all.
3. **Mapping the template name onto the study's copy.** Step 1 yields a CANlab template
   name, but the study's copy is renamed and usually shortened
   (`prep_3a_run_second_level_regression_and_save` becomes
   `cfs_secondlevel_m14_s6_prep_3a_run_regression`), sharing no suffix with it. The two
   are matched on token overlap, requiring the step designator (`prep_3a`, `c2a`, `a2`)
   to agree. Without this a report and the `.mat` files from the same run end up under
   two different names and are documented twice.

In `proj_cfs` that identifies 72 of 98 result files. The remainder are genuinely
ambiguous — `res_warnings.mat` is mentioned by several scripts — and still get their
commits resolved; only the per-script dirty narrowing is unavailable, and they are
marked accordingly.

Intersecting the two recovers the commit that was genuinely checked out. Output is
written as **sidecar files** — existing reports are never modified, since they are
typically already datalad-saved and git-annexed:

```
<model>/results/notes/provenance_<model>.tsv     one row per artifact x dependency
<model>/results/notes/provenance_<model>.mat
<model>/results/html/<script>_provenance.html    ONE PAGE PER RUN (see below)
<model>/results/html/provenance_<model>.html     overview of EVERY artifact in the model
<subdir>/provenance_<proj>_<subdir>.tsv          summary, at the root of the tree it
                                                 documents (secondlevel/, firstlevel/)
```

**One page per run, not per file.** A script usually produces a report *and* one or more
`.mat` files; a page each would say the same thing several times. Artifacts are grouped by
the run that produced them and listed together on one page. The group key is the script
**plus the subject**, where there is one: at second level a script runs once per model, but
at first level the same script runs once per subject — 135 of them in `proj_cfs`, spanning
five different CanlabCore commits — so subject has to be part of the key or every subject
would collapse into a single page. The page is named after the report when the run produced
one, so it sits next to the report; otherwise after the script.

**First-level models have no `results/` directory** — subjects sit directly under the model
with their reports in `sub-*/diagnostics/`. The tooling follows whichever layout is already
there rather than inventing a `results/` tree, so first-level output goes to
`<model>/provenance/` and the per-run pages sit in each `sub-*/diagnostics/`.

### Which evidence is used when a run left several traces

The artifacts in a group are ranked, and the best-evidenced one drives the dependency
table. Strongest first:

| Evidence | Why |
|---|---|
| `mat_header` | `Created on:` inside the `.mat` header. Embedded, carries a time of day, survives git-annex and copying |
| `mat_header+append` | same stamp, but the file was appended to later, so the later write is used |
| `DC.date+mtime` | the report's run *date*, with time of day from the mtime because the two agree |
| `DC.date+eod` | only the date was recoverable — the report is annexed, so its mtime is the annex-add time. End of day is conservative |
| `mtime_only` | no embedded date at all |

The others are listed with their own resolved commits, so when they disagree — which
happens when a later script appends to an existing `.mat` — that is visible rather than
averaged away.

Runs whose script has since been deleted are still resolved, by recovering the source
embedded in the report itself. Those are flagged `script_deleted`, and are exactly the
cases where the record would otherwise be lost entirely.

## Dependency documentation

```matlab
R = LaBGAScore_dep_report('/data/master_github_repos/LaBGAScore', 'index', IDX);
```

Writes `dependencies.tsv`, `dependencies.yml` and `DEPENDENCIES.md` into the
repository root. All three are **generated — never edit them by hand**. The LaBGAS
website collects the `.yml` via `scripts/refresh_dependencies.py`.

## What the fields mean

**Resolution** — how confidently a commit was recovered:

| Value | Meaning |
|---|---|
| `reflog_exact` | run time known to the second — any evidence carrying a time of day (`mat_header`, `mat_header+append`, `DC.date+mtime`) |
| `reflog_day` | only the run *date* was recoverable, because the report is git-annexed and its mtime is the annex-add time. End of day is used, which is conservative |
| `reflog_day_APPROX` | as above, but the repository also moved that same day, so the commit may be one move too late |
| `predates_clone_APPROX` | the run is older than this clone's reflog. Nothing is guessed |
| `no_reflog` | the repository has no reflog at all |

**DirtyStatus** — whether uncommitted edits matter for this script:

| Value | Meaning |
|---|---|
| `clean_for_this_script` | nothing this script reaches is modified today, so the commit hash is a complete description |
| `dirty_now_relevant` | files this script reaches carry uncommitted edits today, named in `RelevantModifiedFiles` |
| `unknown` | no dependency map was available for this artifact — its producing script could not be identified |

**MinDepth** — the shortest path from the script to that repository. `1` means the script
calls it directly; deeper entries are reached through another dependency (CanlabCore calls
MediationToolbox at depth 4, for instance). Sort or filter on this to separate a script's
own dependencies from its dependencies' dependencies.

**Confidence** on a dependency edge — see `DEPENDENCIES.md`, which explains
`resolved` / `ambiguous` / `dotcall` / `dynamic` / `unparseable`.

### Keeping the record free of repositories nothing uses

Six rules stop the call graph inventing dependencies. Each was added after a real false
positive, and together they cut the `proj_cfs` second-level record from 1639 rows to 552,
eliminating BrainSpace, gift, cocoanCORE, MediationToolbox, MKDA, ooFmriDataObjML,
ExploreASL, JuSpace, osprey and `labgas.github.io` entirely. What survives —
CanlabCore, CanlabPrivate, CANlab_help_examples, canlab_single_trials, LaBGAScore,
Neuroimaging_Pattern_Masks, RobustToolbox, CoSMoMVPA — has been traced to real calls:

- **A classdef's property declarations are not calls.** mtree reports them as identifiers,
  so CanlabCore's `@atlas` property `labels` resolved to BrainSpace's `labels.m`.
- **A dot-call may reinforce a repository, never introduce one.** `obj.name` is
  indistinguishable from a property or struct-field read, and properties are the common
  case here (`atlas_obj.labels`, `V.dim`, `V.fname`). You cannot call a method on a class
  you never constructed, and constructing it is a plain call - so if nothing else in the
  file reaches that repository, the dot-call is a field read.
- **An ambiguous name may likewise only reinforce.** Several repositories define
  `spm_defaults` and `get_var`; which one runs is undecidable, so it cannot be evidence
  for a repository the file does not otherwise touch.
- **A core-language builtin is never a third-party call.** Everything under
  `matlabroot/toolbox/matlab` — `isempty`, `find`, `table`, `split`, `error` — is the
  builtin, full stop. CanlabCore defining `@region/table` does not make every `table()`
  call a CanlabCore dependency. Names outside that tree (`predict` lives in `stats`,
  `ttest` in `stats`) are left to the other rules, so `@fmri_data/predict` survives.
- **Evidence for a repository must be repo-unique.** `spm_defaults` exists in spm12,
  cocoanCORE *and* CanlabPrivate; treating it as evidence put all three in play, which is
  how cocoanCORE entered records for scripts that never touch it. Only a name belonging to
  exactly one repository counts.
- **Genuine cross-repository ties are broken with `which()`.** MATLAB's path order is what
  actually decides at run time, so it is asked — but only to pick among candidates the
  index already found, never as the primary resolver. Several copies of a name *inside*
  one repository (spm12 ships two `spm_vol.m`) are labelled `ambiguous_within_repo` and
  are harmless.

Only an unshadowed, repo-unique plain call - `load_atlas`, `merge_atlases`, `fmri_data` -
can introduce a repository. Artifacts whose producing script could not be identified have
no call graph at all, so they fall back to the repositories named in `'followrepos'` rather
than listing every clone on the machine.

## What this cannot tell you

- **Whether a dependency was dirty *at the time*.** Uncommitted work is not timestamped
  anywhere. An edit made before a run and later reverted leaves no trace. What the
  tooling does instead is bound the uncertainty: it names the files this script reaches
  that are modified *today*, which in practice is a short list or empty.
- **Which of several ambiguous candidates actually ran.** Deciding whether `threshold`
  means `@image_vector/threshold` or `@statistic_image/threshold` needs type inference
  that workspace-chained scripts do not support. All candidates are listed.
- **Anything reached through `feval`/`eval`/`str2func`.** Those edges are invisible to
  static analysis; files using them are flagged `dynamic`.

## Design notes

**Nothing shells out.** Commit, branch, remote, reflog and dirty state are all read
from the plain-text files under `.git` (`HEAD`, `refs/`, `packed-refs`, `logs/HEAD`,
`config`, `index`). This is much faster across ~30 repositories, and MATLAB's
`system()`/`!` can fail or hang outright in some sessions, which would otherwise take
a whole analysis down with it. `LaBGAScore_prov_gitstatus` reimplements git's own
dirty-check (compare size and mtime against the index; hash only what differs) and is
verified to reproduce `git status --porcelain` exactly.

**`matlab.codetools.requiredFilesAndProducts` is not used.** It aborts entirely on a
syntax error anywhere in the transitive closure, and there are 33 such files in the
CANlab tree. `LaBGAScore_dep_map` records those as `unparseable` and carries on.

**Resolution is index-based, not `which()`-based.** `which` returns one answer decided
by path order, and on this setup it is the wrong one for the calls that matter most:
`which('predict')` returns a liblinear `.mexa64` from SPM's decoding toolbox, and
`which('ttest')` returns MATLAB's Statistics Toolbox, where the scripts actually call
`@fmri_data/predict` and `@fmri_data/ttest`. The index enumerates every candidate.
`which()` is used only as a *tie-break* among candidates the index already found, where
several repositories genuinely define the same name — there it is exactly the right
question, and it cannot miss a class method because the candidate list is fixed first.

**The transitive closure over-approximates**, because it follows every candidate of an
ambiguous name. For deciding whether an uncommitted edit could have affected a run
that is the safe direction: it flags a file that may not have been used, rather than
missing one that was.

## Files

| File | Purpose |
|---|---|
| `labgascore_prov_protect_reflogs.sh` | disable reflog expiry (run first) |
| `LaBGAScore_prov_protect_reflogs.m` | archive reflogs, report unprotected clones |
| `LaBGAScore_prov_gitinfo.m` | read commit/branch/remote/reflog from `.git` |
| `LaBGAScore_prov_gitstatus.m` | uncommitted-change detection from `.git/index` |
| `LaBGAScore_prov_snapshot.m` | record provenance for a run |
| `LaBGAScore_prov_publish.m` | drop-in `publish` replacement |
| `LaBGAScore_prov_resolve_retrospective.m` | reconstruct provenance for past runs |
| `LaBGAScore_dep_build_index.m` | index every callable file under the dependency roots |
| `LaBGAScore_dep_map.m` | per-script call graph |
| `LaBGAScore_dep_report.m` | generate `DEPENDENCIES.md` / `.tsv` / `.yml` |
| `LaBGAScore_check_display.m` | check whether this session's display can produce good report figures |
