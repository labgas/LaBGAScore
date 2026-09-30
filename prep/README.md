# prep — study setup, BIDS conversion, and pre-first-level preparation

Everything that has to happen between the scanner and
[`firstlevel/`](../firstlevel/README.md): getting images into BIDS, defining the
study's paths, turning stimulus-presentation logfiles into BIDS `events.tsv`
files, and smoothing fMRIPrep output.

Two of the six scripts are generic and run essentially unmodified in any study
(`s0`, `s2`). The other four are **worked examples** that must be read and
rewritten for each study — `s1` in particular, because logfiles differ with
design, task, and stimulus-presentation software. This README says which is
which, and what each script promises the next one.

---

## Contents

- [What this is](#what-this-is)
- [Where prep sits in the workflow](#where-prep-sits-in-the-workflow)
- [Generic vs study-specific](#generic-vs-study-specific)
- [The scripts](#the-scripts)
  - [`LaBGAScore_prep_s0_define_directories.m`](#labgascore_prep_s0_define_directoriesm)
  - [`LaBGAScore_prep_parrec2bids.m`](#labgascore_prep_parrec2bidsm)
  - [`LaBGAScore_prep_s1_write_events_tsv.m`](#labgascore_prep_s1_write_events_tsvm)
  - [`LaBGAScore_prep_s1_write_events_tsv_multisess_multitask.m`](#labgascore_prep_s1_write_events_tsv_multisess_multitaskm)
  - [`LaBGAScore_prep_s2_smooth.m`](#labgascore_prep_s2_smoothm)
  - [`LaBGAScore_prep_s2_smooth_multisess.m`](#labgascore_prep_s2_smooth_multisessm)
- [How the scripts hand off](#how-the-scripts-hand-off)
- [The `events.tsv` contract with firstlevel](#the-eventstsv-contract-with-firstlevel)
- [What lands on disk](#what-lands-on-disk)
- [Dependencies](#dependencies)
- [Known limitations and traps](#known-limitations-and-traps)

---

## What this is

Six MATLAB **scripts** (no functions, no `functions/` subfolder — `prep/` is the
one pipeline domain without one). They are meant to be copied into a study's
`code` subdataset and renamed with that study's prefix, following the convention
in the repo root [`README.md`](../README.md): a study copy of
`LaBGAScore_prep_s0_define_directories.m` becomes
`<study_prefix>_prep_s0_define_directories.m`.

That renaming is not cosmetic. Every script except `parrec2bids` reaches `s0`
through

```matlab
eval([study_prefix '_prep_s0_define_directories']);
```

so the study prefix is what binds a study's scripts together. With
`study_prefix = ''` (the shipped default) that line evaluates
`_prep_s0_define_directories`, which is not a valid MATLAB identifier — the
script stops with *"Unrecognized function or variable"*. **Setting
`study_prefix` is therefore the first edit in every script, not an option.**

## Where prep sits in the workflow

```
scanner
  │
  ├─ PAR/REC ──▶ parrec2bids ─────────────▶ BIDS/       (Philips; study-specific)
  │              (or dcm2bids/heudiconv for other scanners)
  │
  └─ logfiles ─▶ s1_write_events_tsv ─────▶ BIDS/sub-*/func/*_events.tsv
                                             │
BIDS/ ──▶ fMRIPrep (outside this repo) ──▶ derivatives/fmriprep/
                                             │
                                             ▼
                                        s2_smooth ──▶ s6-*.nii.gz
                                             │
                                             ▼
                                    firstlevel/ s1 → s2 → s3
```

`s0` is not a step in that chain — it is called *by* the other steps (and by
`firstlevel`, and by the second-level templates in the
`CANlab_help_examples` LaBGAS fork) to define paths and subject lists.

**Run every script from the root of the superdataset**, e.g.
`/data/proj_discoverie`. `s0` sets `rootdir = pwd`, so starting anywhere else
silently points the whole pipeline at the wrong tree — the single most common
setup mistake in this pipeline.

## Generic vs study-specific

| script | status | what must change per study |
|---|---|---|
| `s0_define_directories` | **generic** | `study_prefix` only (plus `githubrootdir` off the LaBGAS server) |
| `s2_smooth` | **generic** | `study_prefix`; optionally `fwhm`, `prefix`, `subjs2smooth` |
| `s2_smooth_multisess` | **generic** | as above, plus `nr_sess` |
| `parrec2bids` | **example** | scanner protocol, series names, task names, slice timing — 15 `STUDY-SPECIFIC` markers in the code |
| `s1_write_events_tsv` | **example** | essentially all of it: logfile format, column names/types, condition labels, rating parsing |
| `s1_write_events_tsv_multisess_multitask` | **example** | as above, plus session/task structure |

The `s1` scripts are the important case. They read the raw output of whatever
stimulus-presentation software the study used (Presentation, PsychoPy, E-Prime,
…), and **there is no generic version of that problem**. What the shipped scripts
provide is a worked example of the *shape* of the solution — where the loop over
subjects and runs goes, how time zero is established from the first scanner
pulse, how onsets and durations are derived, how the result is written — not code
to be run as-is.

## The scripts

### `LaBGAScore_prep_s0_define_directories.m`

Defines the standard LaBGAS BIDS/DataLad directory layout and the subject lists,
and checks the environment. Called from the superdataset root; called by
everything else.

Defines: `rootdir` (`= pwd`), `githubrootdir` (hardcoded
`/data/master_github_repos`), `sourcedir`, `BIDSdir`, `codedir`, `derivdir`
(`derivatives/fmriprep`), `spmrootdir`, and six subject variables —
`sourcesubjs`/`BIDSsubjs`/`derivsubjs` (names) and
`sourcesubjdirs`/`BIDSsubjdirs`/`derivsubjdirs` (full paths).

Two checks:

1. **SPM12 on the path** — errors with the exact `addpath` command to run if
   `spm.m` is not found. No SPM function is called here; `spmrootdir` is derived
   for later scripts.
2. **The three subject lists match** — errors if `sourcedata`, `BIDS` and
   `derivatives/fmriprep` disagree on which subjects exist. This is deliberately
   strict: a mismatch means preprocessing is incomplete or a subject was added
   without being pushed through the chain, and every later script indexes these
   lists positionally.

It also does `addpath(genpath(codedir),'-end')` if `codedir` is not already on
the path, warning when it does so.

### `LaBGAScore_prep_parrec2bids.m`

Converts Philips PAR/REC to BIDS-organized NIfTI using Xiangruili's
[`dicm2nii`](https://github.com/xiangruili/dicm2nii), and writes the slice-timing
information that fMRIPrep needs for slice-timing correction — derived from the
exam card's slice scan order, fold-over direction and fat-shift direction, which
are script options rather than anything readable from the files.

Study-specific by construction; the shipped version is the `proj_bitter-reward`
two-session, two-task example. It does **not** call `s0` — it runs before a
`derivatives` tree exists, so `s0`'s three-way subject check could not pass.

For non-Philips data, use `dcm2bids` or `heudiconv` instead; nothing in this repo
covers DICOM conversion.

### `LaBGAScore_prep_s1_write_events_tsv.m`

Reads one Presentation `.log` per run from
`sourcedata/sub-*/logfiles/`, extracts onsets, durations and trial-by-trial
ratings, and writes one `events.tsv` per run into `BIDS/sub-*/func/`. Optionally
accumulates all subjects' ratings into a single `phenotype.tsv`.

The shipped example is LaBGAS
[`proj_erythritol_4a`](https://gin.g-node.org/labgas/proj_erythritol_4a): a
four-substance sweet-taste task with delivery, swallow/rinse and rating events.

The parts worth understanding before adapting it:

- **Time zero is the first scanner pulse**, not the logfile's own clock:
  `time_zero = log.Time(log.Trial == 0 & log.EventType == 'Pulse')`, subtracted
  from every timestamp, then divided by 10000 to convert Presentation's
  0.1 ms units to seconds.
- **Durations are differences between successive onsets** — event *n*'s duration
  is `onset(n+1) - onset(n)`. This makes the task's event stream gapless by
  construction, and means the final event of a run cannot get a duration this
  way (see [Known limitations](#known-limitations-and-traps)).
- **Fixation events are dropped** by giving them `NaN` duration and filtering on
  that, so `fixation` acts as the implicit baseline in the first-level model
  rather than as a modelled condition.
- **Ratings are parsed out of the logfile's `Code` strings** by position
  (`scorestring(1,end-3:end)` and friends). This is the most study-specific code
  in the repo and will not survive contact with a different logfile.

### `LaBGAScore_prep_s1_write_events_tsv_multisess_multitask.m`

The same job for designs with more than one session and/or more than one task per
session, writing `sub-*_ses-*_task-*_run-*_events.tsv`. At 1028 lines it is the
largest script in `prep/`, mostly because the per-run body is repeated for each
branch of the subject/session/task structure.

Use this variant when the design has sessions; it is the counterpart of
`firstlevel`'s `s1a`/`s2a` multisession scripts, and the events filenames it
writes are what
`firstlevel/functions/LaBGAScore_firstlevel_find_events.m` resolves by BIDS
inheritance.

### `LaBGAScore_prep_s2_smooth.m`

Unzips fMRIPrep's `*preproc_bold*.nii.gz`, smooths with SPM12, re-zips the
smoothed images, and deletes the unzipped ones. Generic.

Options: `fwhm` (default 6 mm), `prefix` (default `s6-`), `subjs2smooth` (empty
= all subjects). The SPM batch for each subject is saved as
`<sub>_smooth.mat` beside the images, so the smoothing is reproducible from the
batch alone.

`prefix` and `fwhm` must be kept consistent with each other by hand — the prefix
is a label, not something derived from the kernel — and the firstlevel `DSGN`
must glob for the same prefix (`DSGN.funcnames` uses `s6*` in every shipped
example).

### `LaBGAScore_prep_s2_smooth_multisess.m`

The same, looping over `ses-<n>` subdirectories for `nr_sess` sessions. Session
labels are built with `sprintf('ses-%d',sess)`, i.e. **unpadded** (`ses-1`, not
`ses-01`) — datasets using zero-padded session labels need this line changed.

Both smoothing scripts return to `rootdir` when finished, which matters because
`s0` sets `rootdir = pwd`.

## How the scripts hand off

Everything flows through the variables `s0` puts in the base workspace. There is
no struct, no function signature, no return value — the scripts are sequential
and share a workspace, which is why they must be run in order and from the right
directory.

| variable | set by | consumed by |
|---|---|---|
| `rootdir`, `BIDSdir`, `derivdir`, `sourcedir`, `codedir` | `s0` | `s1`, `s2`, all of `firstlevel`, second-level templates |
| `spmrootdir` | `s0` | `firstlevel` `s2` (SPM path checks) |
| `sourcesubjs`, `BIDSsubjs`, `derivsubjs` | `s0` | `s1` (`subjs2write` intersection), `s2` (`subjs2smooth` intersection) |
| `sourcesubjdirs`, `BIDSsubjdirs`, `derivsubjdirs` | `s0` | `s1` (logfile input, events output), `s2` (images), `firstlevel` |
| `study_prefix` | set by hand in each script | the `eval` that calls `s0` |

The three subject-name lists and the three path lists are **positionally
aligned** — `sourcesubjdirs{k}`, `BIDSsubjdirs{k}` and `derivsubjdirs{k}` are the
same participant. `s0`'s three-way equality check is what guarantees that, and
it is why the check is an error rather than a warning.

## The `events.tsv` contract with firstlevel

`firstlevel` reads these files in `s2`/`s2a` and matches columns **by name**, not
position:

- `onset` — seconds from the first scanner pulse of that run
- `duration` — seconds
- `trial_type` — condition label; matched against `DSGN.conditions` with
  `contains`, so a `DSGN` condition name may carry a prefix the events file does
  not (`bit_high calorie` matches `high calorie`). This is how session-specific
  regressor names are built in multisession designs from session-agnostic events
  files.
- any further columns are available as parametric modulators (see
  `LaBGAS_options.pmods` in [`firstlevel/README.md`](../firstlevel/README.md)).
  The `s1` example writes `rating` for this purpose.

Two consequences worth knowing:

- **Column order does not matter** and is not BIDS-canonical in LaBGAS studies
  (BIDS asks for `onset` and `duration` first; the example writes
  `onset, trial_type, rating, duration`, and some studies write `trial_type`
  first). Nothing in the pipeline depends on order, but external BIDS validators
  will comment.
- **Condition names must be `contains`-unambiguous.** If one condition's name is
  a substring of another's, both match the shorter events label and onsets land
  in the wrong regressor with no error. Keep labels distinct, and prefer names
  that differ at the start.

Where the same timings apply to every run, a single inherited file at the BIDS
root (`task-<label>_events.tsv`) is enough —
`LaBGAScore_firstlevel_find_events` searches the run-specific path first, then
walks up to the dataset root. Per-subject files are required whenever timings
vary by subject, which includes any self-paced or jittered design.

## What lands on disk

```
sourcedata/sub-*/logfiles/*.log            input to s1 (never modified)
BIDS/sub-*/func/*_events.tsv               written by s1
BIDS/phenotype/<pheno_name>                written by s1, optional, all subjects
derivatives/fmriprep/sub-*/func/
    *preproc_bold.nii.gz                   fMRIPrep output (input to s2)
    s6-*preproc_bold.nii.gz                written by s2
    <sub>_smooth.mat                       the SPM batch s2 ran
```

`s2` leaves no unzipped `.nii` behind. Under DataLad, run `datalad save` **after**
smoothing finishes, never while a fit or smoothing job is running — saving
annexes the tree and makes files read-only underneath a running job.

## Dependencies

- **SPM12** on the MATLAB path *without* subfolders — checked by `s0`, used by
  `s2` (`spm_select`, `spm_jobman`).
- **LaBGAScore** on the path *with* subfolders.
- **[`dicm2nii`](https://github.com/xiangruili/dicm2nii)** — `parrec2bids` only.
- **fMRIPrep** — run outside MATLAB, between `parrec2bids`/`s1` and `s2`.
- No CANlabCore dependency in `prep/`; that starts at `firstlevel`.

`DEPENDENCIES.md` at the repo root is the generated, authoritative record — do
not hand-edit it (see [`clean/README_provenance.md`](../clean/README_provenance.md)).

## Known limitations and traps

Established by reading the code, not inferred.

### Fixed on 2026-09-30

An audit of these scripts found six real defects, all since repaired. They are
listed because study copies made before that date still carry them:

| script(s) | defect |
|---|---|
| both `s1` | `pheno_tsv = false` stopped on an undefined `pheno_file` — the accumulation and write after the subject loop were not guarded by `if pheno_tsv`, so the documented option did not work |
| both `s1` | the duplicate-logfile check tested `size(logfilepath,1) > 1` on the single path string `fullfile` returns, which can never exceed one row; two matching logfiles fell through to the `~isfile` branch and were reported as *"logfile missing"*, silently skipping the run. The check now runs on the `dir()` result |
| `s1` (single-session) | the "ambiguity about time zero" errors referenced an undefined `subjs{sub}`, so that branch died on the undefined variable instead of reporting the real problem |
| both `s1` | an **empty** `time_zero` was not caught, only a duplicated one; a logfile with no `Trial == 0` `Pulse` row silently yielded an empty table. Now an explicit error |
| both `s2` | `gzip('s6*')` hardcoded the default prefix while `prefix` is an option: changing `prefix` matched nothing, left the smoothed images unzipped, and the following `delete('*.nii')` removed them — **silently destroying the smoothing output**. Both zip and delete are now driven by the list `gunzip` returns |
| `s0` | `spmrootdir` was derived with `strsplit(...,'/spm.m')`, hardcoding a forward slash and so failing on Windows; now `fileparts` |

Three further improvements went in at the same time: the `sourcedata` and `BIDS`
listings are now filtered to directories like the `derivatives` one (a stray file
named `sub-*` used to be counted as a subject); an empty listing is now reported
as *"no sub-\* directories found"* rather than surfacing as a confusing list
mismatch; and `s0`'s `CLEAN UP OBSOLETE VARIABLES` section, previously an empty
heading, now clears what it promises.

The two `s2` scripts were also made consistent with each other: the smoothing
body existed in four verbatim copies (two per file, one per branch of the
`subjs2smooth` test) and is now a single loop in each, over a subject index
resolved beforehand.

### Remaining

**`s1` (both variants)** — the per-run body is still duplicated between the
`subjs2write` and all-subjects branches (and, in the multisession variant,
across the session/task branches), so a fix has to be applied in every copy.
They currently agree; keep them that way. This was left alone deliberately:
unlike the smoothing body, the per-run block differs subtly between branches,
and these are example scripts each study rewrites anyway.

**`s1` (single-session)** — `log.onset(m+1)` and `log.rating(n+3)` are unbounded.
Both encode study-specific facts about this logfile's event ordering: the `+3`
is where the rating sits relative to its event, and the `m+1` works only because
the run is assumed to end on a `fixation` event, which takes the `NaN` branch. A
logfile ending on any other event raises an index error. Left as-is because a
guard would mask a genuine data problem in code every study replaces.

**`s2` (both)** — `prefix` and `fwhm` are independent; keeping them consistent is
still manual, and the firstlevel `DSGN.funcnames` must glob for the same prefix.

**`parrec2bids`** — Philips-only, and study-specific throughout (15
`STUDY-SPECIFIC` markers). Slice timing is reconstructed from exam-card options
that are **not** validated against the data; getting
`slice_scan_order_exam_card` wrong produces a plausible-looking but wrong
slice-timing field that fMRIPrep will apply silently. Only
`fold_over_direction_exam_card == 'AP'` is implemented.

### Static analysis

`clean/LaBGAScore_check_all_scripts.m` (MATLAB's `checkcode`) catches none of the
defects above — they are undefined variables, unreachable branches and ordering
problems, not parse errors. It reports 0 parse-level messages across `prep/`
both before and after these fixes. `clean/use_before_def.py` and
`clean/set_after_use.py` cover the two ordering classes. The rest needs reading,
which is how this list was produced.

### Not yet tested

These fixes are verified by reading and by `checkcode`, **not** by running the
scripts: that needs a study's logfiles and fMRIPrep output. Run `s2` on a single
subject (`subjs2smooth = {'sub-xx'}`) before trusting it on a cohort, and check
that the `s6-*.nii.gz` appear and no stray `.nii` is left behind.
