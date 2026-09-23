# LaBGAScore

Core scripts (and templates for them) for LaBGAS's (Laboratory for Brain-Gut Axis Studies, KU Leuven) standard neuroimaging analysis workflow.

> ### New to the lab? Start here
>
> 1. **[`LaBGAS_fMRI_analysis_workflow.md`](LaBGAS_fMRI_analysis_workflow.md)** — the
>    step-by-step guide to running a study end to end, from server access and DICOM
>    conversion through first and second level. Read it before this README; everything
>    else assumes it. Its **"Ten traps when adapting a template"** section is the
>    single highest-value page in this repo — every trap in it is a real
>    wrong-but-clean run someone already paid for.
> 2. **This README** — what lives where in this repo, and why.
> 3. **[`CANlab_help_examples`](https://github.com/labgas/CANlab_help_examples)** (LaBGAS
>    fork) — the second-level templates you will actually copy into your study. See
>    [Relationship to `CANlab_help_examples`](#relationship-to-canlab_help_examples-labgas-fork).
>
> Never edit a checked-in template in place for a study. Copy it into your project's
> `code` subdataset and adapt the copy — see
> [Repository structure and naming convention](#repository-structure-and-naming-convention).

## Contents

- [What this is](#what-this-is)
- [Typical LaBGAS project structure](#typical-labgas-project-structure)
- [Repository structure and naming convention](#repository-structure-and-naming-convention)
- [Dependencies](#dependencies)
- [Domain-by-domain overview](#domain-by-domain-overview)
- [Relationship to `CANlab_help_examples` (LaBGAS fork)](#relationship-to-canlab_help_examples-labgas-fork)
- [Running scripts and publishing reports](#running-scripts-and-publishing-reports)
- [Provenance and dependency documentation](#provenance-and-dependency-documentation)
- [Tests, CI, and the script checkers](#tests-ci-and-the-script-checkers)
- [License](#license)

## What this is

LaBGAScore is a curated collection of MATLAB scripts and helper functions implementing LaBGAS's standard neuroimaging analysis workflow, spanning BIDS conversion, first-level and second-level fMRI modeling, MVPA/machine-learning pipelines, PET, MRS, and several auxiliary toolkits. A small amount of code is in other languages — SAS macros under `stats_tools/sas_macros/`, and standalone Python utilities under `qr_code/`. Scripts can either be run directly from this repo, or — more commonly for a real study — copied into that study's own project repo and adapted there for study-specific purposes.

This is not a packaged software product: there is no build system or automated test suite. Three lightweight checkers stand in for a test suite (see [Tests, CI, and the script checkers](#tests-ci-and-the-script-checkers)).

## Typical LaBGAS project structure

Every LaBGAS project/study is organized as a **DataLad superdataset** `proj_xxx`, with subdatasets:

- `sourcedata` — raw/source acquisitions
- `BIDS` — BIDS-converted imaging and phenotype data
- `derivatives` — preprocessing output (e.g. `derivatives/fmriprep/`)
- `code` — analysis scripts, adapted per study
- `firstlevel` — first-level model outputs
- `secondlevel` — second-level (group) model outputs

LaBGAScore's own scripts are canonical/example copies: for a real study they are downloaded and adapted into that study's `code` subdataset, not run in place from this GitHub repo. `firstlevel/LaBGAScore_firstlevel_s1_options_dsgn_struct.m`'s own header states this explicitly ("provided in LaBGAScore as an example and needs to be downloaded and adapted to the code subdataset for your study/project"), citing the real study repo `proj_erythritol_4a` as a worked example.

For second-level (group) analyses specifically, there is a matching pair of folders per model: `code/secondlevel/model_x/` (the adapted scripts) and `secondlevel/model_x/` (that model's outputs — `masks/`, `results/`, `results/figures/`, `results/notes/`, `results/html/`). This pairing is set up by the sibling `CANlab_help_examples` repo's second-level template scripts (see [below](#relationship-to-canlab_help_examples-labgas-fork)), not by anything in LaBGAScore itself.

## Repository structure and naming convention

The top level is organized by analysis topic, not as a conventional MATLAB toolbox:

`prep/`, `firstlevel/`, `secondlevel/`, `stats_tools/`, `atlas_mask_tools/`, `pet/`, `mrs/`, `decoding_toolbox/`, `cosmomvpa/`, `graphvar/`, `juspace/`, `power/`, `figures/`, `clean/`, `qr_code/`.

Most domains split into `<domain>/scripts/` (or top-level `.m` scripts) and `<domain>/functions/` (helper functions the scripts call). There is exactly one class in the whole repo, `secondlevel/classes/ProgressTracker.m` (a `handle`-derived progress/ETA helper for long loops).

**Naming conventions:**
- Numbered/lettered prefixes (`prep_1_`, `prep_2_`, `s0_`, `s1_`, `s2_`, `s3_`, `a_`, `a2_`, and lettered variants like `s1a_`/`s1b_`/`s2a_`/`s2b_`) mark ordered, reusable pipeline steps meant to be copied into a study's `code` subdataset and adapted — not generic library code.
- A `_example` suffix (e.g. `secondlevel/scripts/LaBGAScore_secondlevel_ooFmriDataObjML_example.m`) marks a script that is illustrative only, not meant to be copied and run as-is.
- Scripts (files with no `function`/`classdef` declaration) use a standard MATLAB comment header: `%% scriptname.m` title, a `*USAGE*` section, optionally `*OPTIONS*`/`*DEPENDENCIES*`/`*NOTES*`, then an author/date/version block. Functions follow MATLAB's standard function help-text convention (H1 line, syntax) instead.

**Study model directories.** Inside a study (not in this repo), a second-level model is named
after the **first-level model it is built on**:

```
model_<firstlevel number><letter>_<description>
```

so `model_2a_casecontrol_cov_scanner` is built on `firstlevel/model_2_basic`, is the first
second-level analysis derived from it, and controls for scanner. Scripts inside follow
`<proj>_secondlevel_m<N><letter>_s<step>_<template name>.m`. The first-level model is
otherwise invisible from the second-level name, and it is the thing most likely to change
underneath an analysis, so two results are only comparable if they share it. Models named
before this convention keep their names — renaming would break paths recorded in saved
`.mat` files and published reports. Full detail, including how to tell what an
older model was built on, is in
[`LaBGAS_fMRI_analysis_workflow.md`](LaBGAS_fMRI_analysis_workflow.md#naming-convention).

`secondlevel/` additionally contains seven detailed, standalone usage guides — `README_ENet_neuroimaging_pipeline.md`, `README_ENet_plotting.md`, `README_PLSDA_neuroimaging_pipeline.md`, `README_PLSDA_paired_neuroimaging_pipeline.md`, `README_PLSDA_plotting.md`, `README_PLSR_neuroimaging_pipeline.md`, `README_PLSR_plotting.md` — for the Elastic Net / PLS-DA / PLSR pipeline functions and their diagnostic plotting companions. Consult those directly for usage details; they aren't duplicated here.

## Dependencies

Most workflows require **[CanlabCore](https://github.com/canlab/CanlabCore)**, **[Neuroimaging_Pattern_Masks](https://github.com/canlab/Neuroimaging_Pattern_Masks)**, and **[SPM12](https://www.fil.ion.ucl.ac.uk/spm/)** on the MATLAB path. None of these is vendored into this repo — all must be cloned/installed separately and added to the path (see `prep/LaBGAScore_prep_s0_define_directories.m` and `pet/scripts/LaBGAScore_pet_a_set_up_paths_always_run_first.m`, which check for and set up these paths per study).

Several domain folders assume their own external toolbox, also not vendored, expected pre-installed and on the path:

| Folder | External toolbox |
|---|---|
| `cosmomvpa/` | [CoSMoMVPA](https://www.cosmomvpa.org/) |
| `decoding_toolbox/` | [The Decoding Toolbox (TDT)](https://sites.google.com/site/tdtdecodingtoolbox/) |
| `graphvar/` | [GraphVar](http://rfmri.org/GraphVar) |
| `juspace/` | [JuSpace](https://github.com/juryxy/JuSpace) |
| `mrs/` | [Osprey](https://github.com/schorschinho/osprey) |
| `stats_tools/sas_macros/` | [SAS](https://www.sas.com/) with SAS/STAT (not MATLAB — these are `.sas` macro files) |

Exception: `pet/functions/LCN_*.m` are legacy KU Leuven PET-processing functions vendored directly into this repo rather than pulled from an external toolbox.

There is no package manager or dependency manifest — path setup is manual, via the scripts below.

**Path/environment setup scripts:**
- `prep/LaBGAScore_prep_s0_define_directories.m` — checks SPM12 is on the path (errors with setup instructions if not) and does `addpath(genpath(codedir),'-end')` for a study's `code` directory.
- `pet/scripts/LaBGAScore_pet_a_set_up_paths_always_run_first.m` — the PET-domain equivalent; must be run first before other LaBGAScore PET second-level batch scripts.
- `clean/LaBGAScore_move_repos_matlabpath.m` — reorders two repos on the MATLAB path and calls `savepath`.
- `clean/LaBGAScore_smart_parallel_pool_setup.m` — configures a `parpool`, capping worker count at ~60% of available cores.

Both assume local repos live under `/data/master_github_repos` (see `githubrootdir` in `prep_s0`).

## Domain-by-domain overview

**`prep/`** — BIDS conversion and pre-first-level prep: `LaBGAScore_prep_parrec2bids.m` (Philips PAR/REC → BIDS-organized NIfTI), `LaBGAScore_prep_s0_define_directories.m` (defines study paths, checks SPM12), `LaBGAScore_prep_s1_write_events_tsv*.m` (logfiles → BIDS `events.tsv`, single- or multi-session/task), `LaBGAScore_prep_s2_smooth*.m` (unzip → smooth → re-zip fMRIPrep output).

**`firstlevel/`** — SPM/CANlab first-level GLM pipeline: `s1_options_dsgn_struct` (builds the CANlab-style `DSGN` design struct) → `s2_fit_model` (fits and diagnoses first-level models, cloning CANlab dependencies as needed) → `s3_diagnose_model` (publishes an HTML diagnostic report via CANlab's `scn_spm_design_check`, saves VIFs), with `s1a`/`s2a` as the multisession/multitask variant of the first two. A separate phMRI (pharmacological challenge) chain (`s1b`/`s2b`/`s3b`/`s3c`) models the post-administration period as one regressor per timebin and builds SPM batches directly rather than through `canlab_glm_*`. `functions/canlab_glm_subject_levels*_old.m` are older CANlab GLM functions kept because current example scripts still call them. See [`firstlevel/README.md`](firstlevel/README.md) for the task-fMRI chain — the `LaBGAS_options` reference, the noise model, the inter-script handoff, what lands on disk, and the contrast-order contract with the second-level scripts — and [`firstlevel/README_phMRI.md`](firstlevel/README_phMRI.md) for the phMRI chain.

**`secondlevel/`** — group-level statistics and MVPA/ML pipelines: TFCE permutation inference using classic TFCE, Smith & Nichols 2009 (`group_tfce_from_subject_maps.m` → `tfce_one_fmri_dat.m` → `tfce_volume.m` → `tfce_transform_3d.m`), PLS-DA/PLSR/Elastic Net pipelines with matching diagnostic-plotting functions (see the seven `README_*.md` guides above), the shared per-fold helpers they are built from (`foldPreprocess.m`, `residualizeFold.m`, `residualizeY.m`, `applyScaling.m`, `capLV.m`, `validateCovariates.m`, `warnUnknownOptions.m`, `setParforStream.m`, `globalBaselineCV.m`, `selectENetHyperparams.m`, `enetLambdaGrid.m`, `logitSafe.m`, `quickCV_*.m`, `bootstrapOOB_*.m`, `makeGroupedFolds.m`, `swapWithinSubjectLabels.m`, `quickGroupedCV.m`), atlas/threshold validation helpers (`validateAtlasLabels.m`, `maskToSignificant.m`), dice-overlap tools (`dice_statistic_image*.m`), blob-reporting helpers that give any thresholded map the montage/region-table treatment `c2a` gives a GLM result (`LaBGAScore_blob_montage.m`, and `LaBGAScore_pdm_report.m` / `LaBGAScore_pdm_report_image.m` for PDM mediation results, with `LaBGAScore_pdm_regenerate_reports.m` to backfill that reporting from already-written `PDM*.nii` without refitting the bootstrap), an ROI/parcel extraction script (`LaBGAScore_secondlevel_extractparcels_sessions.m`), a MACS-toolbox model-space batch-setup script (`LaBGAScore_secondlevel_MS_mat_pipeline.m`), an MVPA-regression-on-connectivity-betas script (`LaBGAScore_secondlevel_mvpa_beta_maps_conn.m`), a PLS/ENet ROI-pipeline wrapper (`LaBGAScore_secondlevel_roi_run_plot_PLS_ENet_pipeline.m`), an object-oriented ML toolkit example (`LaBGAScore_secondlevel_ooFmriDataObjML_example.m`), and the `ProgressTracker` class.

**`stats_tools/`** — general-purpose statistics helpers, in two languages. `functions/LaBGAScore_Storey_FDR.m` implements several multiple-comparison corrections behind one interface, selected with `'method'`: Storey q-values (`'sas'`, the default, reproducing SAS PROC MULTTEST's PFDR — spline then bootstrap on SAS's own trigger; also `'lambda'`, `'spline'`), plain Benjamini-Hochberg (`'bh'`), the two adaptive FDR procedures (`'adaptivefdr'` = Benjamini & Hochberg 2000 with the lowest-slope pi0, `'bky'` = Benjamini, Krieger & Yekutieli 2006 two-stage), and FWER step-down corrections (`'stepdown_sidak'`, `'holm'`). A reliability guard can reject an implausible pi0 and fall back to BH; it is ON for the pi0-estimating methods but OFF for `'sas'`, because SAS has no such check and a q-value that PROC MULTTEST cannot reproduce defeats the purpose of a SAS mode. Override with `'guard'`. The function's header carries the measured behaviour of each estimator, the comparison against R's `qvalue` package, and what the original Storey papers do and do not say about the number of tests. `LaBGAScore_combat_fit.m` and `LaBGAScore_combat_apply.m` split ComBat into a fit step and an apply step, so harmonisation parameters can be estimated on a training fold only and applied to held-out data — the form ComBat has to take inside cross-validation, where fitting on all the data would leak. `LaBGAScore_dummy_code.m` expands a phenotype column into k-1 indicator columns, which is what an unordered factor such as scanning site requires in a design matrix. `sas_macros/` holds SAS macros for statistics the MATLAB side does not cover: `mixed_effectsize.sas` (`%mixed_effectsize`) computes effect sizes — partial eta-squared with a noncentral-F confidence interval, Cohen's f², and optionally eta-squared/omega-squared — for the fixed effects of a model fitted with `PROC MIXED`, and `es_identify.sas` (`%es_identify`, `%es_identify_ds`) audits effect sizes in an existing results table to determine which statistic was actually reported. Both files are self-contained; `%include` the one you need. They require SAS with SAS/STAT, and are the only non-MATLAB code in the analysis workflow proper. See [`stats_tools/sas_macros/README.md`](stats_tools/sas_macros/README.md) for usage, the formulas used, the marginal-versus-conditional residual distinction, and the caveats to state in a Methods section — including its note that the macros have not yet been run against a real model in SAS.

**`atlas_mask_tools/`** — `LaBGAScore_atlas_binary_mask_from_atlas.m` and `LaBGAScore_atlas_rois_from_atlas.m` generate custom atlas/mask and per-ROI objects from a chosen atlas; the folder also ships a set of ready-made brain/gray-matter mask and template NIfTIs (`brain_masks/`, `brain_templates/`, `gray_matter_masks/`).

**`pet/`** — PET-specific pipeline: path setup (`a_set_up_paths_always_run_first.m`, `a2_set_default_options.m`), DICOM-to-BIDS conversion, preprocessing, kinetic modeling (e.g. `LaBGAScore_pet_model_TSPO_DPA714.m`), prep/signature-application scripts mirroring the fMRI `prep_*` naming, and a PLS/ENet pipeline runner. `functions/LCN_*.m` are the vendored legacy KU Leuven PET-processing functions.

**`mrs/`** — MRS pipeline built on Osprey: BIDS conversion (`LaBGAScore_prep_mrs2bids.m`), GE/Philips Osprey jobfiles (single- and multi-session), and run scripts (`LaBGAScore_mrs_run_osprey_GE.m`/`_Philips.m`).

**`cosmomvpa/`** — `LaBGAScore_cosmomvpa_searchlight_rsa.m`: first-level representational similarity analysis (behavioral vs. neural dissimilarity, leave-one-run-out cross-validation) plus group-level permutation/TFCE testing, via CoSMoMVPA.

**`decoding_toolbox/`** — TDT-based decoding: `LaBGAScore_decoding_template_xclass_acc.m` (first-level (cross-)classification accuracy plus group-level testing) and `LaBGAScore_decoding_SVM_between_subjects.m` (a full between-subject SVM decoding pipeline with permutation-based TFCE inference, sharing `secondlevel/functions/tfce_volume.m` with the second-level TFCE stack, with optional in-fold ComBat harmonisation and nuisance residualisation). `LaBGAScore_export_scaled_contrasts.m` writes second-level contrast objects out as one NIfTI per subject, so a decoding analysis can be run on exactly the images the GLM used rather than on raw first-level contrasts.

**`graphvar/`** — `LaBGAScore_prep_graphvar_input_from_conn.m`: builds GraphVar input files from CONN toolbox ROI-to-ROI connectivity output.

**`juspace/`** — `LaBGAScore_prep_juspace_input.m` (builds first-level con-image input for JuSpace PET-receptor spatial correlation analysis) and `LaBGAScore_juspace_corr_behav.m` (correlates behavioral data with JuSpace's spatial-correlation results).

**`power/`** — standalone power-analysis helper functions: `holmthreshold.m`, `medianIQR_to_meanSD.m`, `print_LEAR_matrix.m`, `sd_from_ci.m`.

**`figures/`** — plotting helpers: `canlabCmap.m` (CANlab colormap), `cluster_surface_plots.m`, `save_all_open_figures_smart.m`.

**`clean/`** — utilities: `LaBGAScore_clean_gzip_all_nii.m` (recursive `.nii` gzip), `LaBGAScore_clean_sourcedata.m` (cleans the `sourcedata` subdataset), `LaBGAScore_move_repos_matlabpath.m`, `LaBGAScore_smart_parallel_pool_setup.m`, the provenance/dependency tooling and the report runners described below (`labgascore_run_headless.sh`, `LaBGAScore_run_reports.m`), and the three script checkers — `LaBGAScore_check_all_scripts.m` plus the two Python option-ordering checkers `use_before_def.py` and `set_after_use.py`, with their positive controls in `checker_positive_controls/` (see [Tests, CI, and the script checkers](#tests-ci-and-the-script-checkers)).

**`qr_code/`** — standalone Python utilities (`emailQR.py`, `email_QR_Input.py`, `QRtoPDF.py`) for generating and emailing QR codes; not MATLAB, and unrelated to the neuroimaging pipeline proper.

## Relationship to `CANlab_help_examples` (LaBGAS fork)

[`CANlab_help_examples`](https://github.com/labgas/CANlab_help_examples) is a sibling repo whose `Second_level_analysis_template_scripts/` directory provides LaBGAS's group-level (second-level) fMRI analysis templates — univariate GLM, cross-validated SVM, CANlab signature-pattern responses, searchlight correlation, and single-trial/runwise MVPA and mediation. Several of those scripts depend directly on LaBGAScore:

| LaBGAScore function/script | Called from (in CANlab_help_examples) | Purpose |
|---|---|---|
| `LaBGAScore_prep_s0_define_directories` (`prep/`) | `a_set_up_paths_always_run_first.m` | Defines `rootdir`, `codedir`, `BIDSdir`, `spmrootdir`, etc. for the study |
| `LaBGAScore_firstlevel_s1_options_dsgn_struct` (`firstlevel/`) | `a_set_up_paths_always_run_first.m` | Builds the first-level `DSGN` struct the second-level framework reads `DSGN.modeldir`/`DSGN.conditions`/`DSGN.contrastnames` from |
| `LaBGAScore_firstlevel_s2_fit_model.m` (`firstlevel/`) | upstream of `prep_3f_create_fmri_data_single_trial_object.m` (not called directly) | Fits first-level models; produces the single-trial con images that single-trial/runwise MVPA and mediation analyses consume |
| `LaBGAScore_atlas_binary_mask_from_atlas.m` (`atlas_mask_tools/`) | referenced via `atlasname_glm`/`atlasname_svm` options in `a2_set_default_options.m` | Generates custom `.mat` atlas/mask objects usable as GLM/SVM masks |
| `LaBGAScore_atlas_rois_from_atlas.m` (`atlas_mask_tools/`) | referenced via `roi_names`/`roi_modelname`/`roi_set_name` options | Generates per-ROI atlas objects for ROI-average analysis |

In short: LaBGAScore owns study setup, first-level modeling, and atlas/mask generation; `CANlab_help_examples` (LaBGAS fork) owns the second-level/group analysis templates built on top of LaBGAScore's outputs. For that repo's own internals, see its own [`README.md`](https://github.com/labgas/CANlab_help_examples/blob/master/Second_level_analysis_template_scripts/README.md) under `Second_level_analysis_template_scripts/`.

## Running scripts and publishing reports

**Headless is the default.** Analysis scripts are run from the Linux command line with no
display, using the wrapper in `clean/`:

```bash
labgascore_run_headless.sh -d /data/proj_xxx \
    -s proj_secondlevel_m1_s0_a_set_up_paths_always_run_first \
    -a data_objects.mat \
    proj_secondlevel_m1_s4_prep_2_load_image_data_and_save
```

Prepend `setsid nohup ... > run.log 2>&1 < /dev/null &` for long chains, which then survive
logout. `-h` prints the full option list. Run interactively in X2go instead only when you
want higher-resolution figures (72 dpi headless against 96–144 dpi in a graphical session)
or are debugging. Figure *size* is not limited headless: without a display `publish` prints
figures rather than capturing them from the screen.

Two traps worth knowing if you write the MATLAB invocation yourself:

- **`matlab -batch` cannot `publish`** — it fails with *"Unable to run the `publish`
  function, because it is not supported for this ..."*. Use `-nodisplay` with `-r`.
- **Redirect stdin from `/dev/null`**, or MATLAB may exit before running anything.

**`clean/LaBGAScore_run_reports.m`** does the publishing and, more importantly, decides
whether it worked. `publish` catches a script's error into the html report and returns
normally, so a crashed run raises no exception and exits 0 — a chain of scripts can appear
to complete while one of them died, sometimes after an hour of computation and before
anything was saved. `LaBGAScore_run_reports` reads each report back, fails on the error
markup `publish` writes (`<pre class="codeoutput error">`, matched as markup rather than by
searching the report text, since phrases like "Error in" occur in ordinary comments), and
optionally asserts that the results file the script should have written exists and is of
plausible size. It continues past failures and names them all at the end.

```matlab
results = LaBGAScore_run_reports({'script_one' 'script_two'}, htmlsavedir, ...
              'artefacts', {'data_objects.mat' []});
if ~all(results.ok), error('%d report(s) failed', sum(~results.ok)); end
```

Full instructions, including the X2go DPI table for the interactive route, are in
[`LaBGAS_fMRI_analysis_workflow.md`](LaBGAS_fMRI_analysis_workflow.md).

## Provenance and dependency documentation

Scripts are copied into a study's `code` subdataset and frozen there, but their
**dependencies are not**: CanlabCore, our fork of CANlab_help_examples, CanlabPrivate and
the rest are shared clones that keep moving. `clean/` holds tooling that closes that gap.
See **[`clean/README_provenance.md`](clean/README_provenance.md)** for the full guide.

**Record provenance when you publish a script** — replace the `publish` call documented in
each script header:

```matlab
LaBGAScore_prov_publish('my_secondlevel_script', htmlsavedir)
```

The html report gains a Provenance section listing the commit of every dependency the
script reaches, plus the screen and figure dimensions it was produced at, and a
machine-readable copy lands in the model's `results/notes/`.

**Reconstruct provenance for analyses already run** with
`LaBGAScore_prov_resolve_retrospective`, which recovers the same information from each
artifact's own embedded date and each clone's git reflog. It writes sidecar files and
never modifies the existing ones.

It covers **result `.mat` files as well as published reports** — only a few scripts per
model are ever published, so reports alone would leave most of the pipeline undocumented.
Output is one page per *run*, grouping a script's report with the `.mat` files it wrote and
ranking them by how well each is dated.

**Check your display before publishing.** `publish()` captures figures from the screen, so
your X2go window size and DPI decide how figures in the report come out:

```matlab
LaBGAScore_check_display        % what can this session produce, and what to change
```

Recommended settings per screen size are in
[`LaBGAS_fMRI_analysis_workflow.md`](LaBGAS_fMRI_analysis_workflow.md#2-set-up-your-x2go-display-for-publishing-figures).

> **Run [`clean/labgascore_prov_protect_reflogs.sh`](clean/labgascore_prov_protect_reflogs.sh)
> once per machine.** Git prunes reflogs after 90 days by default, silently. The reflog is
> the only record of which commit a clone actually had checked out at a past moment, and it
> cannot be reconstructed afterwards.

**[`DEPENDENCIES.md`](DEPENDENCIES.md)** documents what every script in this repo calls and
which repository each of those lives in. It, `dependencies.tsv` and `dependencies.yml` are
**generated** by `clean/LaBGAScore_dep_report.m` — regenerate them rather than editing.

## Tests, CI, and the script checkers

There is no automated test suite and no CI pipeline in this repo. Verify changes manually by running the affected script against real or study data; several scripts support MATLAB's `publish()` to generate a date-stamped HTML report of their output (e.g. `cosmomvpa/LaBGAScore_cosmomvpa_searchlight_rsa.m`, `decoding_toolbox/LaBGAScore_decoding_template_xclass_acc.m`).

What exists instead is three cheap checkers in `clean/`, each aimed at a failure
mode the previous one cannot see. Run all three before launching a long chain —
together they take seconds, and each was written after a specific failure that
cost hours.

| checker | catches | run it on |
|---|---|---|
| `LaBGAScore_check_all_scripts.m` | parse errors, plus style/performance suggestions listed separately | this repo, or any subtree you pass as `rootdir` |
| `clean/use_before_def.py` | an option **read above the line that defines it** | a study's model script directory |
| `clean/set_after_use.py` | an option **set below the line that already consumed it** | the same |

**`LaBGAScore_check_all_scripts.m`** wraps MATLAB's built-in Code Analyzer
(`checkcode`). It reliably catches parse errors — a file that cannot run past the
flagged line — but does **not** catch undefined variables used at runtime, calls
to functions that don't exist or aren't on the path, or logic bugs; those still
require reading the code.

**`use_before_def.py`** covers the gap immediately next to that one. A guarded
default (`if ~exist('opt','var'), opt = ...; end`) only protects an option if it
runs *before* every read. Put the guard beside one use and miss an earlier one,
and the script dies on *"Unrecognized function or variable"* at the earlier line
— which may sit near the end of a long script, after all the expensive work and
before anything is saved. `checkcode` cannot see this: the line is syntactically
perfect. This is not hypothetical — `contrast_objects_tag` was guarded beside
`savefilenamedata` while a `printhdr` eight lines above also read it, and `prep_3`
died there after 75 minutes of ComBat and contrast formation, having saved
nothing.

**`set_after_use.py`** covers the opposite, and worse, case: the option *is*
defined before use, so nothing errors — but the author's real setting sits
*below* the consumer, so the value in force is the default and the setting
silently does nothing. The motivating case: a decoding template builds its output
directory from `results_tag` a few lines after the guard, and a study copy that
set `results_tag` further down had eleven runs write into the same untagged
folder and overwrite each other's maps. The statistics were unaffected, because
results are always recomputed; the saved artefacts were not.

`set_after_use.py` is **advisory, not a gate**. One benign pattern still trips it
— a variable legitimately reused for successive outputs (`savefilename` for a
second file, `figtitle` for the next figure). Expect a handful per model, read
them rather than chase them, and keep the checker because the failure it does
catch is silent, expensive and otherwise invisible.

Both Python checkers take one argument, the directory of scripts to scan:

```bash
python3 clean/use_before_def.py /data/proj_xxx/code/secondlevel/model_N_name
python3 clean/set_after_use.py  /data/proj_xxx/code/secondlevel/model_N_name
```

Neither needs MATLAB. Both ship with a **positive control** — a deliberately
broken miniature script each checker must flag — so a silent PASS can be
trusted rather than assumed:

```bash
clean/checker_positive_controls/run_controls.sh
```

Each control is positive for its own checker and negative for the other, and the
runner asserts both directions, so a checker that starts flagging everything
fails the controls just as a checker that stops flagging anything does. The
checkers exit non-zero when they find something, so in a control run a non-zero
exit is success; the runner translates that and prints a single verdict line.

`checkcode` reports **zero** messages on either control file. That is the point:
both failure modes are invisible to MATLAB's own static analysis, which is why
these two checkers exist alongside it rather than inside it.

See also **["Ten traps when adapting a template"](LaBGAS_fMRI_analysis_workflow.md)**
in the workflow document, which is the catalogue these checkers were distilled
from; traps 6 and 7 are exactly what the two Python scripts automate.

## License

GNU General Public License v3.0 — see [`LICENSE`](LICENSE).
