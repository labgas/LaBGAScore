# firstlevel — LaBGAS first-level fMRI analysis

Reference documentation for the task-fMRI first-level scripts in this folder. It documents
what LaBGAScore itself adds — the `LaBGAS_options` structure, the noise model, the
inter-script handoff, and what lands on disk — and links out to CANlab for the `DSGN`
structure and the SPM machinery underneath, rather than restating them.

The **pharmacological-challenge (phMRI) chain** (`s1b`/`s2b`/`s3b`/`s3c`) is a separate
family that does not use `canlab_glm_*` at all; it has its own guide,
[`README_phMRI.md`](README_phMRI.md).

For the **step-by-step procedure** — creating the `firstlevel` DataLad subdataset,
downloading and renaming these scripts into your study's `code` subdataset, running them,
and saving the output with `datalad save` — see the *First-level analysis* section of
[`LaBGAS_fMRI_analysis_workflow.md`](../LaBGAS_fMRI_analysis_workflow.md). That procedure
is not repeated here.

## Contents

- [What this is](#what-this-is)
- [The scripts](#the-scripts)
- [How the scripts hand off](#how-the-scripts-hand-off)
- [`LaBGAS_options`](#labgas_options)
- [The `DSGN` structure](#the-dsgn-structure)
- [Contrast order is a cross-repo contract](#contrast-order-is-a-cross-repo-contract)
- [The noise model](#the-noise-model)
- [Parametric modulators](#parametric-modulators)
- [Single-trial models](#single-trial-models)
- [What lands on disk](#what-lands-on-disk)
- [Provenance](#provenance)
- [The multisession/multitask variant](#the-multisessionmultitask-variant)
- [`functions/`](#functions)
- [Dependencies](#dependencies)
- [Known limitations](#known-limitations)

## What this is

These scripts fit subject-level (first-level) GLMs on fMRIPrep-preprocessed, smoothed BOLD
data using SPM12 through CANlab's `canlab_glm_*` batch tools, and produce the contrast
images that group-level analysis consumes.

Where they sit in the LaBGAS workflow:

```
prep/                          firstlevel/  (this folder)         CANlab_help_examples
                                                                  (LaBGAS fork)
prep_s0_define_directories  →  s1_options_dsgn_struct         →   Second_level_analysis_
prep_s1_write_events_tsv    →  s2_fit_model                       template_scripts/
prep_s2_smooth              →  s3_diagnose_model
```

Upstream, `prep/` converts to BIDS, writes `events.tsv` files, and smooths the fMRIPrep
output. Downstream, the second-level templates live in the LaBGAS fork of
[`CANlab_help_examples`](https://github.com/labgas/CANlab_help_examples/tree/master/Second_level_analysis_template_scripts)
and read this folder's `DSGN` structure and `con_*.nii` images directly — see
[Contrast order is a cross-repo contract](#contrast-order-is-a-cross-repo-contract).

Like everything in LaBGAScore, these are **templates to copy into a study's `code`
subdataset and adapt**, not library code. Editing a script here changes the template for
future studies, not any running analysis.

## The scripts

| Script | Role | Study-specific? |
|---|---|---|
| [`LaBGAScore_firstlevel_s1_options_dsgn_struct.m`](LaBGAScore_firstlevel_s1_options_dsgn_struct.m) | Sets `LaBGAS_options` and builds the CANlab `DSGN` structure: conditions, parametric modulators, contrasts, timing | **Yes** — heavily. The version in the repo is a worked example from [`proj_erythritol_4a`](https://gin.g-node.org/labgas/proj_erythritol_4a) |
| [`LaBGAScore_firstlevel_s2_fit_model.m`](LaBGAScore_firstlevel_s2_fit_model.m) | Does the work: dependency checks, directory structure, noise regressors, onsets/durations/pmods, design plots, model estimation, contrasts, diagnostics | No — only the two generic script names at the marked lines need renaming |
| [`LaBGAScore_firstlevel_s3_diagnose_model.m`](LaBGAScore_firstlevel_s3_diagnose_model.m) | Publishes the per-subject diagnostic HTML report (design, collinearity, VIFs, thresholded contrast montages) | No — and **not for standalone use**; `s2` publishes it |
| [`LaBGAScore_firstlevel_s1a_options_dsgn_multisess_multitask.m`](LaBGAScore_firstlevel_s1a_options_dsgn_multisess_multitask.m) | `s1` for designs with more than one session and/or more than one task | **Yes** |
| [`LaBGAScore_firstlevel_s2a_fit_model_multisess_multitask.m`](LaBGAScore_firstlevel_s2a_fit_model_multisess_multitask.m) | `s2` for the same | No |

`s1`/`s2` and `s1a`/`s2a` are alternatives, not steps: pick the pair that matches your
design. Both pairs publish the same `s3`.

## How the scripts hand off

These are **scripts, not functions**. They take no arguments and return nothing; all state
travels through the base MATLAB workspace, and each step chains to the previous one through
an `exist` guard:

- `s2` checks `if ~exist('DSGN','var')` and runs `s1` itself when the design structure is
  not already in the workspace
  ([`s2:132`](LaBGAScore_firstlevel_s2_fit_model.m#L132)).
- `s1` (and `s1a`, `s1b`) call
  [`prep/LaBGAScore_prep_s0_define_directories.m`](../prep/LaBGAScore_prep_s0_define_directories.m),
  which derives `rootdir`, `BIDSdir`, `derivdir`, `codedir`, `spmrootdir` and the subject
  lists from the current working directory.
- `s2` publishes `s3` at the end of each subject's iteration, which is why `s3` can rely on
  `subjfirstdir`, `DSGN` and `LaBGAS_options` already existing.

Three consequences worth internalising:

1. **You run `s2` only.** It runs `s1` and `s3` for you.
2. **Run it with the superdataset root as the working directory** (e.g.
   `/data/proj_erythritol/proj_erythritol_4a`). `prep_s0` derives every path from `pwd`;
   run from anywhere else and the directory detection is wrong.
3. **Variable names in `prep_s0` are an interface.** Renaming `rootdir`, `derivdir`,
   `BIDSdir` or the subject-list variables silently breaks every script downstream of it,
   in this repo and in the second-level fork.

## `LaBGAS_options`

LaBGAScore's own options structure, set at the top of `s1`/`s1a` and consumed by `s2`/`s3`.
The values below are the LaBGAS defaults as shipped; they may be study-specific, so discuss
with Lukas if in doubt.

### `LaBGAS_options.mandatory`

| Field | Default | Meaning |
|---|---|---|
| `spike_def` | `'fMRIprep'` | Source of spike regressors. `'fMRIprep'` uses fMRIPrep's `motion_outlier` columns (thresholds set by fMRIPrep's `--fd-spike-threshold`/`--dvars-spike-threshold`). `'CANlab'` uses [`scn_session_spike_id`](https://github.com/canlab/CanlabCore/blob/master/CanlabCore/diagnostics/scn_session_spike_id.m) (Mahalanobis distance) plus DVARS |
| `omit_spike_trials` | `'no'` | Whether to drop trials whose onset coincides with a spike. **Not recommended** — LaBGAS prefers to handle this later, from VIFs in single-trial models |
| `spikes_percent_threshold` | `0.15` | Proportion of volumes flagged as spikes above which a warning is printed |
| `vif_thresh` | `2` | Variance inflation factor above which `s3` flags a regressor. Sensible range 1.3 (stringent) to 5 (lenient) |

`movement_reg_quadratic` (default `true`) sits directly on `LaBGAS_options`, not under
`mandatory`; set it `false` to omit the quadratic terms of the motion parameters and their
derivatives.

### `LaBGAS_options.subjs2analyze` (optional)

Cell array of subject directory names in `derivdir` to analyse, e.g.
`sort({'sub-01','sub-02'})`. Leave it empty (or comment it out) to loop over every subject.
A name that is not present in `derivdir` raises an error naming the offending subject.

### `LaBGAS_options.spikes`

| Field | Default | Meaning |
|---|---|---|
| `dvars_threshold` | `2` | Standardized DVARS above which a volume counts as a spike. **Required** when `spike_def = 'CANlab'`; `s1` errors if it is missing. CANlab's own default is 3, which LaBGAS considers lenient |
| `spike_additional_vols` | `0` | Number of volumes after each spike to additionally regress out. **Not recommended**: the approach comes from resting-state work, and it creates missingness that is not at random in task data |

### `LaBGAS_options.pmods`

Required only if your design has parametric modulators. See
[Parametric modulators](#parametric-modulators).

| Field | Default | Meaning |
|---|---|---|
| `pmod_polynom` | `1` | Polynomial expansion of the modulator; 2 or 3 for quadratic/cubic (mind orthogonalization) |
| `pmod_name` | `'rating'` | Column name of the modulator in the `events.tsv` files |
| `pmod_ortho_off` | `false` | `true` demeans modulators per condition, bypassing SPM's serial orthogonalization. Set it when you have more than one modulator per condition |
| `pmod_type` | `'parametric_standard'` | `'parametric_standard'` gives each condition an unmodulated *and* a modulated regressor; `'parametric_singleregressor'` gives only the modulated one. Both are options of [`onsets2fmridesign`](https://github.com/canlab/CanlabCore/blob/master/CanlabCore/Model_building_tools/onsets2fmridesign.m) |

### `LaBGAS_options.display`

Plotting, and thresholding/masking of the first-level maps **for display only** — these
never affect the fitted model.

| Field | Default | Meaning |
|---|---|---|
| `plotdesign` | `true` | Save per-run design plots. Turning it off saves time; not recommended |
| `plotmontages` | `true` | Save thresholded montages of the contrast maps in the report. Turning it off saves a lot of time; not recommended |
| `input_threshold` | `0.005` | p-value or raw-value range, depending on `thresh_type` |
| `thresh_type` | `'unc'` | One of `'fdr'`, `'bfr'`, `'unc'`, `'extent'`/`'cluster_extent'`, `'raw-between'`, `'raw-outside'` — passed straight to CANlab's `statistic_image.threshold` |
| `k` | `25` | Cluster extent threshold, in voxels |
| `mask` | `which('gm_mask_canlab2023_coarse_fmriprep20_0_20.nii')` | Image to mask the displayed maps with |

## The `DSGN` structure

`DSGN` is **CANlab's** structure, not LaBGAS's, and it is documented properly upstream.
Rather than duplicating that here, use these:

| Source | How to read it | What it gives |
|---|---|---|
| [`canlab_glm_dsgninfo.txt`](https://github.com/canlab/CanlabCore/blob/master/CanlabCore/GLM_Batch_tools/canlab_glm_dsgninfo.txt) | `canlab_glm_subject_levels('dsgninfo')` | **The authoritative field-by-field reference** — every required and optional field, its default, and an example |
| [`canlab_glm_README.txt`](https://github.com/canlab/CanlabCore/blob/master/CanlabCore/GLM_Batch_tools/canlab_glm_README.txt) | `canlab_glm_subject_levels('README')` | Overview of the `canlab_glm_*` family and what must exist before you run it |
| [`canlab_glm_example_DSGN_setup.txt`](https://github.com/canlab/CanlabCore/blob/master/CanlabCore/GLM_Batch_tools/canlab_glm_example_DSGN_setup.txt) | — | Worked example, including the flexible contrast-specification syntax |
| [`promptDSGN.m`](https://github.com/canlab/CanlabCore/blob/master/CanlabCore/GLM_Batch_tools/promptDSGN.m) | `promptDSGN` | Michael Sun's interactive `DSGN` builder — the easiest way to start a new model |
| [`canlab_spm_contrast_job_luka.m`](https://github.com/canlab/CanlabPrivate) (CanlabPrivate) | `help canlab_spm_contrast_job_luka` | Contrast specification in detail |

Further background: CANlab's [first-level SPM12
tutorial](https://canlab.github.io/_pages/tutorials/html/first_level_spm12.html), and the
`CANlabReposGuide_Hackpad.pdf` in the LaBGAS [Google
Drive](https://drive.google.com/drive/folders/1-M5UvibmsWXVCIrR31-qJNu506pDA_0t).

### Where LaBGAS deviates from CANlab/SPM defaults

| Field | LaBGAS value | Why |
|---|---|---|
| `hpf` | `180` | SPM's default is 128 s. CANlab uses 180 s because responses to sustained stimuli (pain, food) last long, and variance is lost at shorter filter lengths |
| `ar1` | left off (SPM default is `true`) | Tor Wager's recommendation: the AR(1) algorithm pools across the whole brain and behaves poorly in some situations. Less concerning when the endpoint is a group analysis |
| `fmri_t` / `fmri_t0` | number of slices / middle slice | Set for multiband acquisition **with** slice timing correction. SPM's and CANlab's defaults (16 / 1) are appropriate for multiband *without* slice timing |
| `notimemod` | `true` | No linear time modulation of conditions |
| `regmatching` | `'regexp'` | Contrasts are matched by regular expression against the beta names in `SPM.Vbeta.descrip` — see below |

### The contrast-matching idiom

With `DSGN.regmatching = 'regexp'`, each contrast is written as a regular expression
matched against beta regressor names. The idiom used throughout the shipped example:

```matlab
DSGN.contrasts{c} = {{'.*sucrose{1}\s[^x]'}};   % unmodulated regressor only
DSGN.contrasts{c} = {{'.*liking_sucrose'}};     % modulated regressor only
```

`\s[^x]` matches a single whitespace *not* followed by `x`, which selects the unmodulated
regressor and excludes the `sucrose x liking_sucrose` interaction term SPM generates for
the parametric modulator. A two-element cell array gives a difference contrast, weighted
`[1 -1]`.

## Contrast order is a cross-repo contract

The second-level templates in the `CANlab_help_examples` LaBGAS fork consume this folder's
output directly, and the coupling is tighter than it looks:

- `a_set_up_paths_always_run_first.m` runs your `..._firstlevel_s1_options_dsgn_struct`
  script if `DSGN` is not in the workspace, then takes
  `[~,modelname] = fileparts(DSGN.modeldir)` and sets
  `datadir = fullfile(rootdir,'firstlevel',modelname)`. **The model directory name is the
  join key between the two repos.**
- `prep_1_set_conditions_contrasts_colors.m` sets `DAT.conditions` from
  `DSGN.contrastnames`, and locates images as
  `fullfile(datadir, DAT.subfolders{i}, DAT.functional_wildcard{i})` — in practice
  `firstlevel/<model>/sub-*/con_0001.nii`, `con_0002.nii`, …

> **The order of `DSGN.contrasts` fixes the `con_000N` numbering that second level
> hardcodes.** Inserting or reordering a contrast in `s1` silently remaps every downstream
> second-level analysis onto the wrong images — no error, just different results. If you
> change the contrast list, either append to the end, or re-check
> `DAT.functional_wildcard` and `DAT.conditions` in the second-level scripts.

## The noise model

`s2` writes one `noise_regs.mat` per run, containing a matrix `R` in the format SPM and the
CANlab GLM tools expect. By default it holds:

1. **24 head motion parameters** — the six realignment parameters, their first derivatives,
   and the quadratic terms of both, standardized. The derivatives are recomputed from the
   parameters rather than taken from fMRIPrep, and the quadratics are computed last so they
   stay orthogonal. Set `movement_reg_quadratic = false` for the 12-column version.
   Rationale: CANlab's [first-level SPM12
   tutorial](https://canlab.github.io/_pages/tutorials/html/first_level_spm12.html).
2. **The global CSF signal**, taken straight from fMRIPrep's `csf` confound column.
3. **Spike regressors**, per `spike_def` — either fMRIPrep's `motion_outlier` columns, or
   CANlab's Mahalanobis-distance spikes plus non-redundant DVARS spikes. Redundant spike
   columns (a volume flagged by more than one criterion) are collapsed to one.

The proportion of volumes flagged as spikes is reported per run, with a warning above
`spikes_percent_threshold`.

## Parametric modulators

Supply `DSGN.pmods` (one cell per session, one cell per condition) and set
`LaBGAS_options.pmods`. `s2` reads the modulator column named by `pmod_name` from each
`events.tsv`, and computes three versions of it:

- **raw**;
- **demeaned per run** — useful in some designs, and unlike SPM's behaviour;
- **demeaned per condition** — what SPM does internally.

Which one enters the model follows from `pmod_type` and `pmod_ortho_off`. SPM serially
orthogonalizes multiple modulators within a condition, which is rarely what you want with
more than one modulator per condition; setting `pmod_ortho_off = true` sets `orth = 0` and
uses the per-condition demeaned values instead. Background:
[Bob Spunt on parametric modulation](https://www.bobspunt.com/resources/teaching/single-subject-analysis/parametric-modulation/),
the [MRC CBU wiki](https://imaging.mrc-cbu.cam.ac.uk/imaging/ParametricModulations), and
`help onsets2fmridesign`.

Design plots for the modulated design are written alongside the unmodulated ones.

## Single-trial models

Setting `DSGN.singletrials` (a cell per session, a logical per condition) or
`DSGN.singletrialsall = true` converts the flagged conditions into one regressor per trial.
`s2` then builds per-trial contrasts by calling
`canlab_spm_contrast_job_single_trials_lukasvo` — which lives in **CanlabPrivate**, not
CanlabCore.

These per-trial `con_*` images are what
`prep_3f_create_fmri_data_single_trial_object.m` consumes in the second-level fork, feeding
the single-trial/runwise MVPA and mediation analyses.

## What lands on disk

Two trees. In the **derivatives** subdataset, per subject and run:

```
derivatives/fmriprep/sub-xx/func/
├── run-1/
│   ├── s6-*_bold.nii                     unzipped smoothed image, copied here by s2
│   ├── *_desc-confounds_timeseries.tsv   fMRIPrep confounds, copied here by s2
│   └── <modelingfilesdir>/
│       ├── noise_regs.mat                the R matrix described above
│       ├── <condition>.mat               onsets, durations, pmods per condition
│       ├── design_sub-xx_run-1.png       design matrix plot
│       └── design_<pmod_type>_sub-xx_run-1.png
└── run-2/ …
```

In the **firstlevel** subdataset, per subject:

```
firstlevel/<model>/
├── sub-xx/
│   ├── SPM.mat
│   ├── beta_*.nii, con_*.nii, spmT_*.nii     con_000N order = DSGN.contrasts order
│   ├── vifs.mat                              variance inflation factors, saved by s3
│   └── diagnostics/
│       ├── LaBGAScore_firstlevel_s3_diagnose_model.html   the published report
│       └── *.png
└── provenance/
    └── sub-xx/                               one snapshot per subject, see below
        ├── *_provenance_<timestamp>.tsv
        └── *_provenance_<timestamp>.mat
```

`s2` creates all of these directories with `mkdir` if they do not exist — including
`firstlevel/` itself. **Create `firstlevel/` as a DataLad subdataset first**, as the
workflow document describes: a plain `mkdir` gives you an ordinary directory inside the
superdataset rather than a subdataset, and the imaging output then lands in the wrong place
in the DataLad hierarchy.

## Provenance

Scripts are frozen per study in the `code` subdataset, but their dependencies are not:
CanlabCore, CanlabPrivate and SPM12 are shared clones that keep moving. Which commit of
each produced a given first-level model is therefore not recoverable from the model itself.

`s2`/`s2a` close that gap by publishing the diagnostic report through
[`clean/LaBGAScore_prov_publish.m`](../clean/LaBGAScore_prov_publish.m) rather than a bare
`publish()`. You get the same report, plus:

- a **Provenance section in the HTML**, listing the commit of every dependency the run
  reached, which of them carried uncommitted local changes relevant to this script, and the
  MATLAB and SPM versions in play;
- a machine-readable `.tsv`/`.mat` copy under
  `firstlevel/<model>/provenance/<subject>/`;
- the **screen geometry, DPI and resulting figure dimensions**, with a flag on any figure
  whose size was set by the display rather than by the script.

**One snapshot per subject per model.** Unlike second level, where a script runs once per
analysis, `s3` runs once per subject — hence the per-subject subdirectory. The location
mirrors what
[`clean/LaBGAScore_prov_resolve_retrospective.m`](../clean/LaBGAScore_prov_resolve_retrospective.m)
already uses at first level: `<modeldir>/provenance/`, not the `results/notes` tree it uses
at second level, because first-level models have no `results/` directory.

`s2`/`s2a` pass `maxHeight 800` / `maxWidth 1600`, keeping report figures the size they have
always been here. `LaBGAScore_prov_publish` does not set them by default, because `publish()`
applies them by resizing the PNG on disk — detail is lost, not merely displayed smaller.
Drop those two arguments at the call site for full-resolution montages at the cost of larger
files.

This requires **LaBGAScore itself on the MATLAB path with subfolders**, since
`prov_publish` lives in `clean/`. Full detail in
[`clean/README_provenance.md`](../clean/README_provenance.md).

## The multisession/multitask variant

Use `s1a`/`s2a` when subjects were scanned in more than one session, or performed more than
one task. What changes:

| | `s1`/`s2` | `s1a`/`s2a` |
|---|---|---|
| Sessions | one | `nr_sess`; `s2a` nests a `for ses = 1:nr_sess` loop inside the subject loop, and expects `ses-01`, `ses-02`, … under each subject |
| Tasks | one, implicit | `tasknames` lists them; `taskname` selects **one**. Image, confound and events discovery is filtered by it |
| Run dir names | `rundirnames = {'run-1',…}` | append the task, e.g. `{'run-1_FID','run-2_FID'}` |
| `DSGN.funcnames` | one entry per run | one entry per session × run, with the session and task in the wildcard |

**One task per run of the script, one model per task.** To analyse a second task, change
`taskname`, point `DSGN.modeldir` and `DSGN.modelingfilesdir` at a different model
directory, and run again.

## `functions/`

Two vendored copies of older CANlab functions:

- `canlab_glm_subject_levels_old.m`
- `canlab_glm_subject_levels_run1subject_old.m`

`s2`/`s2a` call the `_old` version deliberately. The current
[`canlab_glm_subject_levels.m`](https://github.com/canlab/CanlabCore/blob/master/CanlabCore/GLM_Batch_tools/canlab_glm_subject_levels.m)
changed how missing onset and duration files are handled, and errors out on the way these
scripts write them ([`s2:942`](LaBGAScore_firstlevel_s2_fit_model.m#L942)). Keep the call
pointed at the vendored copy unless you have re-tested against the current upstream one.

## Dependencies

| Dependency | Needed for | Note |
|---|---|---|
| [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) | everything | On the MATLAB path **without** subdirectories |
| [CanlabCore](https://github.com/canlab/CanlabCore) | `canlab_glm_*`, `onsets2fmridesign`, `scn_spm_design_check`, `scn_session_spike_id`, `statistic_image`, `fmri_mask_image`, `canlab_results_fmridisplay` | `s2` clones it into `githubrootdir` if absent |
| [CanlabPrivate](https://github.com/canlab/CanlabPrivate) | `canlab_spm_contrast_job_luka`, `canlab_spm_contrast_job_single_trials_lukasvo` | Access-restricted. `s2` clones it if absent, which fails without access — needed for contrast generation |
| [LaBGAScore](https://github.com/labgas/LaBGAScore) itself, **with subfolders** | `prep/LaBGAScore_prep_s0_define_directories.m`, and `clean/LaBGAScore_prov_publish.m` for the diagnostic report | Not auto-added; `addpath(genpath(...))` it yourself |
| fMRIPrep output | confounds, preprocessed BOLD | Produced upstream of this repo |
| `events.tsv` files | onsets, durations, modulators | Written by `prep/LaBGAScore_prep_s1_write_events_tsv*.m` |

`s2` shells out to `git clone` for the two CANlab repos when it does not find them under
`githubrootdir` and adds them to the path with `addpath(genpath(...))`. It is not a
side-effect-free script.

## Known limitations

- **`githubrootdir` is hardcoded** to `/data/master_github_repos` in `s1`/`s1a` and
  `prep_s0`. Change it if your clones live elsewhere.
- **Report figures are captured from the screen.** `publish()` grabs figures as they are
  drawn, so your session geometry and display DPI decide how they come out.
  `clean/LaBGAScore_check_display.m` tells you whether the current X2go session can produce
  a full-size figure, and the recommended settings per screen are under *Set up your X2go
  display for publishing figures* in
  [`LaBGAS_fMRI_analysis_workflow.md`](../LaBGAS_fMRI_analysis_workflow.md). Every report
  records the screen geometry, DPI and resulting figure dimensions — see
  [Provenance](#provenance).
- **No tests.** Static analysis over this folder only, from the MATLAB prompt:
  ```
  addpath(genpath('..')); LaBGAScore_check_all_scripts(pwd)
  ```
  `checkcode` catches parse errors but not undefined variables, missing functions, or wrong
  indexing — the defect classes that actually occur here. Read the code path.
- **Two defects were fixed on 2026/09/02** that affected whether a model ran at all:
  designs with parametric modulators errored in the single-session path, and under
  `LaBGAS_options.subjs2analyze` the `s2a` subject loop either fit nothing or errored.
