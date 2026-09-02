# firstlevel — phMRI (pharmacological challenge) analysis

Reference documentation for the pharmacological-fMRI chain in this folder:
`s1b`/`s2b`/`s3b`/`s3c`. These scripts are a **separate family** from the task-fMRI chain
documented in [`README.md`](README.md) — read that one first for the shared context
(directory setup, the `prep_s0` handoff, the noise model, how to run these scripts from the
superdataset root), and this one for everything specific to phMRI.

The most important structural difference: **this chain does not use CANlab's
`canlab_glm_*` tools or the `DSGN` structure at all.** It builds SPM batches directly
through `spm_jobman`. The CANlab `DSGN` documentation linked from `README.md` does not
apply here; `phDSGN` is LaBGAS's own structure and is documented in full below.

## Contents

- [What phMRI modelling means here](#what-phmri-modelling-means-here)
- [The scripts](#the-scripts)
- [The `phDSGN` structure](#the-phdsgn-structure)
- [The `sessions` field is spelled three different ways](#the-sessions-field-is-spelled-three-different-ways)
- [`s1b` — model specification, estimation, contrasts](#s1b--model-specification-estimation-contrasts)
- [`s2b` — covariate contrasts by spline interpolation](#s2b--covariate-contrasts-by-spline-interpolation)
- [`s3b` / `s3c` — signature responses and neurotransmitter maps](#s3b--s3c--signature-responses-and-neurotransmitter-maps)
- [Dependencies](#dependencies)
- [Known limitations](#known-limitations)

## What phMRI modelling means here

In a pharmacological challenge there is no task and there are no events. A substance is
administered part-way through a long continuous scan, and the question is how the BOLD
signal evolves over the minutes that follow.

The design these scripts build reflects that directly: **the post-administration period is
cut into fixed-length timebins, and each timebin becomes one condition regressor.**

```
volume 1                     start_dyn                                    nr_dyns
   |............................|=======|=======|=======| ... |=======|      |
        baseline (not modelled)   bin 1   bin 2   bin 3        bin N
                                  <-- timebin_length seconds each -->
```

From `s1b`:

```matlab
phDSGN.timebin_length_tr = phDSGN.timebin_length / phDSGN.tr;
phDSGN.nr_timebins       = (phDSGN.nr_dyns - phDSGN.start_dyn) / phDSGN.timebin_length_tr;
```

Each bin is specified in **scans**, not seconds (`timing.units = 'scans'`), with onset
`timebin_length_tr * b + (start_dyn - timebin_length_tr)` and duration
`timebin_length_tr`. One session per condition (e.g. sucrose / erythritol / water), so the
model has `nr_timebins × nr_sessions` regressors of interest, named
`<condition>_bin_<b>`.

Two consequences follow, and both differ sharply from task fMRI:

- **The high-pass filter must be very long** — 3540 s in the shipped example, against 180 s
  for task designs. A pharmacological response is slow, and a standard filter would remove
  it as drift.
- **Spikes are not modelled.** The noise model is the 24 head motion parameters plus the
  global CSF signal, and nothing else (`nr_noise_reg = 25` in `s2b`). `phDSGN` has no
  `spike_def` field.

`s1b` sets `cvi = 'AR(1)'` (SPM's default) and `mthresh = -inf`, i.e. no implicit masking.

## The scripts

| Script | Role |
|---|---|
| [`LaBGAScore_firstlevel_s1b_fit_phMRI_model.m`](LaBGAScore_firstlevel_s1b_fit_phMRI_model.m) | The whole first level in one script: noise regressors, model specification, estimation, and between-condition contrasts per timebin. Supports two or three conditions; more would need hardcoded changes |
| [`LaBGAScore_firstlevel_s2b_phMRI_covar_contrasts.m`](LaBGAScore_firstlevel_s2b_phMRI_covar_contrasts.m) | Adds a second set of contrasts to the estimated model, weighted by spline-interpolated covariates (e.g. blood hormones) |
| [`LaBGAScore_firstlevel_s3b_phMRI_signature_responses.m`](LaBGAScore_firstlevel_s3b_phMRI_signature_responses.m) | Applies CANlab signatures to every first-level beta and writes a long-format table for mixed-model analysis |
| [`LaBGAScore_firstlevel_s3c_phMRI_neurotransmitter_maps.m`](LaBGAScore_firstlevel_s3c_phMRI_neurotransmitter_maps.m) | Same, using cosine similarity with Hansen et al. (2022) PET neurotransmitter maps |

All four are **study-specific** and need adapting. `s2b`, `s3b` and `s3c` each require
`s1b` to have run first; `s3b` and `s3c` are independent of each other and of `s2b`.

Unlike the task-fMRI chain, these scripts do **not** chain to each other automatically —
run them in order yourself.

## The `phDSGN` structure

Declared at the top of `s1b`. Every field is study-specific unless noted.

### Not phMRI-specific

| Field | Example | Meaning |
|---|---|---|
| `modelingfilesdir` | `'model_1_basic'` | Subfolder of `firstlevel/` where model output goes, and the name of the model-specific subdirectory created under each session's `func` dir |
| `tr` | `2.5` | Repetition time, seconds |
| `t` | `46` | Microtime resolution — set to the number of slices if slice timing was performed |
| `t0` | `23` | Microtime onset — the reference slice used in slice timing correction |
| `hpf` | `3540` | High-pass filter length, seconds. Very long for phMRI (see above) |
| `multireg` | `'noise_regs'` | Name of the `.mat` file holding the noise regressor matrix, written per session |

### phMRI-specific

| Field | Example | Meaning |
|---|---|---|
| `nr_dyns` | `1416` | Number of dynamics (volumes) per run |
| `start_dyn` | `264` | Volume number at which administration starts |
| `timebin_length` | `60` | Length of each post-administration timebin, seconds |
| `sessions` | `{'sucrose','erythritol','water'}` | Condition names, **in the order of the corresponding sessions** |

### Derived — do not edit

| Field | Computed as |
|---|---|
| `timebin_length_tr` | `timebin_length / tr` |
| `nr_timebins` | `(nr_dyns - start_dyn) / timebin_length_tr` |

`s1b` also fills `phDSGN.conditions` and `phDSGN.subjects` at run time, and saves the whole
structure as `DSGN.mat` next to each subject's `SPM.mat` — which is how `s2b`, `s3b` and
`s3c` recover it.

`s1b` additionally uses a subset of `LaBGAS_options` — `mandatory.omit_spike_trials`,
`mandatory.spikes_percent_threshold`, `mandatory.vif_thresh`, `movement_reg_quadratic`,
`spikes.dvars_threshold`, `spikes.spike_additional_vols` — with the meanings documented in
[`README.md`](README.md#labgas_options).

## The `sessions` field is spelled three different ways

`phDSGN` is **re-declared by hand** in `s1b` and `s2b`, and again as loose variables
(`study_prefix`, `modelingfilesdir`, `sessions`) in `s3b` and `s3c`. Each script's header
says to copy it from `s1b`. Nothing checks that the copies agree, and one field in
particular must be spelled differently in different scripts:

| Script | `sessions` must match | Shipped example |
|---|---|---|
| `s1b` | the condition names you want in the model | `{'sucrose','erythritol','water'}` |
| `s2b` | the `Condition` column of the **covariate `.csv`**, capitalization included | `{'Sucrose','Erythritol','Water'}` |
| `s3b` / `s3c` | the `SPM.Vbeta.descrip` field written by `s1b` | `{'sucrose','erythritol','water'}` |

`s3b` lowercases the covariate file's `Condition` column on read (`covars_table.Condition =
lower(...)`), which is why its `sessions` is lowercase while `s2b`'s is not.

**How it fails:** a mismatch produces no error. `covars_table.Condition == sessions{s}`
simply matches nothing, the covariate vector comes back empty, and the affected contrasts
are silently skipped (`s2b`) or filled with `NaN` (`s3b`/`s3c`). Check the number of
contrasts and the proportion of `NaN` in the output table before trusting a run.

If you change `modelingfilesdir`, `nr_dyns`, `start_dyn`, `timebin_length` or `tr` in one
script, change them in all four — `nr_timebins` is recomputed independently in each, and a
disagreement silently misaligns the contrast indexing.

## `s1b` — model specification, estimation, contrasts

One `spm_jobman` run per subject, in three batches: specification, estimation, contrasts.

**Input discovery.** For each session `ses-<n>` under the subject's derivatives directory,
it expects exactly one smoothed image (`s6*.nii.gz`, gunzipped in place) and exactly one
fMRIPrep confounds file (`*desc-confounds_timeseries.tsv`). More than one of either is an
error; none is a warning, and that session is treated as missing and dropped — so a subject
with a missing condition still gets a model, with fewer sessions.

**Contrasts.** Pairwise between-condition t-contrasts, per timebin: with three sessions it
writes 3 × `nr_timebins` contrasts (session 1 vs 2, 1 vs 3, 2 vs 3), named
`<cond_a>_<cond_b>_bin_<a>`. A subject left with a single session is skipped with a warning
after estimation, so the model exists but has no contrasts. The contrast code is written
for **at most three sessions**.

**Output**, per subject, under `firstlevel/<modelingfilesdir>/sub-xx/`: `SPM.mat`,
`DSGN.mat` (the `phDSGN` structure), `beta_*.nii`, `con_*.nii`, `spmT_*.nii`. The noise
regressor `.mat` files go next to the imaging data, under
`derivatives/fmriprep/sub-xx/ses-<n>/func/<modelingfilesdir>/`.

## `s2b` — covariate contrasts by spline interpolation

Adds contrasts to the already-estimated model in which the **weights are a physiological
or behavioural covariate**, interpolated onto the timebin grid — so the contrast asks where
BOLD tracks the covariate's time course rather than where two conditions differ.

**Input:** a `.csv` in `BIDSdir`, named by `covars_filename`, with columns:

| Column | Purpose |
|---|---|
| `ID` | Subject identifier; matched against the last three characters of the subject directory name |
| `Condition` | Session/condition name — see [the spelling table](#the-sessions-field-is-spelled-three-different-ways) |
| `Timepoint` | Sampling time. A value of `-10` is recoded to `0` before interpolation |
| `Scan_included` | Logical; rows that are false are dropped |
| one per covariate | Named in `name_covars`, e.g. `delta_GLP1` |

**Interpolation:** `interp1(..., 'spline')` onto `1:nr_timebins`, per subject and session,
and only when that session has **more than two** non-NaN values.

**Inclusion rule:** a subject enters the analysis only if it has **at least six** non-NaN
values for each covariate, on the assumption that covariates are either all present or all
missing. Subjects failing this are dropped from `subs2include`, and the `covars` struct is
filtered to match — so the two stay index-aligned.

**Other options:** `nr_noise_reg` (25 in the shipped example: 24 head motion parameters plus
CSF) tells the contrast builder how many nuisance columns to skip when placing weights. Get
this wrong and the weights land on the wrong regressors, with no error.

The script shells out to `git annex unannex *.mat` before writing, because the `SPM.mat`
it modifies is already annexed.

## `s3b` / `s3c` — signature responses and neurotransmitter maps

Both take the first-level betas (one per timebin per condition), reduce each to a single
number per brain map, and write a **long-format table for mixed-model analysis** in the
statistics package of your choice.

| | `s3b` | `s3c` |
|---|---|---|
| Maps | CANlab signatures, from `path_sigs` — either a `.nii` on the path via `which(...)`, or a `load_image_set` keyword | Hansen et al. (2022) PET neurotransmitter maps, `load_image_set('hansen22')`, filtered to `name_nts` |
| Metric | `apply_mask(..., 'pattern_expression')` | cosine similarity |
| Naming | `name_sigs` — list every signature separately, expanding any `load_image_set` set into its constituent names | derived from the map metadata (`target` + `primary_reference`) |

**Output:** `secondlevel/<modelingfilesdir>/SignatureResults_<results_suffix>.xlsx`, one row
per beta, with columns `PPID`, `beta_number`, `beta_descrip`, `condition`, `timebin`, then
one column per map, then — if `add_covars` is true — one column per covariate, interpolated
exactly as in `s2b`.

Note that these write into the **secondlevel** subdataset, since the table is group-level
input. Create that subdataset with DataLad first; the scripts error if it does not exist.

For the mixed models these tables feed, and in particular for reporting effect sizes from
`PROC MIXED`, see [`stats_tools/sas_macros/README.md`](../stats_tools/sas_macros/README.md).

## Dependencies

| Dependency | Needed by | Note |
|---|---|---|
| SPM12 | all | `spm_jobman`, `cfg_dep`, `spm_select` |
| CanlabCore | `s3b`, `s3c` | `load_image_set`, `apply_mask`, `fmri_data`, `fmri_data_st` |
| MATLAB Report Generator | `s2b` | `mlreportgen.utils.capitalizeFirstChar` |
| git-annex on the system path | `s2b` | Shells out to `git annex unannex` |
| MATLAB Statistics Toolbox | `s2b`, `s3b`, `s3c` | `interp1` spline interpolation is base MATLAB; `readtable`/`writetable` are base |

`s3b` and `s3c` locate directories by calling your study's own `prep_s0` script indirectly:

```matlab
eval([study_prefix '_prep_s0_define_directories']);
```

so `study_prefix` must match the renamed copy in your `code` subdataset (e.g. `'ery_ph'` →
`ery_ph_prep_s0_define_directories.m`). `s1b` and `s2b` call
`LaBGAScore_prep_s0_define_directories` by name instead, and that line needs renaming by
hand.

## Known limitations

- **At most three conditions.** `s1b`'s contrast section and `s2b`'s contrast indexing are
  written for two or three sessions; more requires hardcoded changes.
- **`phDSGN` is copied by hand across four scripts** with no consistency check. See
  [the spelling table](#the-sessions-field-is-spelled-three-different-ways).
- **Silent skipping.** Missing sessions, missing covariate values, and condition-name
  mismatches all produce warnings or nothing at all rather than errors. Check the number of
  contrasts written per subject, and the `NaN` count in the results table.
- **Not verified against a live SPM run.** The 2026/09/02 revision fixed several defects in
  this chain (below) by reading and by targeted MATLAB tests, not by running the pipeline
  end to end on real data. Confirm on the LaBGAS server against a real phMRI model before
  relying on it.
- **Results predating 2026/09/02 differ.** That revision fixed: `s2b`'s contrast loop,
  which started at subject index 3 and so silently skipped the first two included
  subjects; `s3c`'s results table,
  which errored outright unless exactly five neurotransmitter maps were selected; a
  `secondlevel` directory check that tested the `firstlevel` directory instead, so `mkdir`
  created `secondlevel/` outside DataLad; a subject count taken from the derivatives list
  rather than the first-level list; and an `fprintf` typo that turned a "no noise files"
  warning into an error.
