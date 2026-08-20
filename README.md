# LaBGAScore

Core scripts (and templates for them) for LaBGAS's (Laboratory for Brain-Gut Axis Studies, KU Leuven) standard neuroimaging analysis workflow.

## Contents

- [What this is](#what-this-is)
- [Typical LaBGAS project structure](#typical-labgas-project-structure)
- [Repository structure and naming convention](#repository-structure-and-naming-convention)
- [Dependencies](#dependencies)
- [Domain-by-domain overview](#domain-by-domain-overview)
- [Relationship to `CANlab_help_examples` (LaBGAS fork)](#relationship-to-canlab_help_examples-labgas-fork)
- [Tests and CI](#tests-and-ci)
- [License](#license)

## What this is

LaBGAScore is a curated collection of MATLAB scripts and helper functions implementing LaBGAS's standard neuroimaging analysis workflow, spanning BIDS conversion, first-level and second-level fMRI modeling, MVPA/machine-learning pipelines, PET, MRS, and several auxiliary toolkits. Scripts can either be run directly from this repo, or — more commonly for a real study — copied into that study's own project repo and adapted there for study-specific purposes.

This is not a packaged software product: there is no build system or automated test suite. A lightweight static-analysis helper is available (see [Tests and CI](#tests-and-ci)).

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

`secondlevel/` additionally contains six detailed, standalone usage guides — `README_ENet_neuroimaging_pipeline.md`, `README_ENet_plotting.md`, `README_PLSDA_neuroimaging_pipeline.md`, `README_PLSDA_plotting.md`, `README_PLSR_neuroimaging_pipeline.md`, `README_PLSR_plotting.md` — for the Elastic Net / PLS-DA / PLSR pipeline functions and their diagnostic plotting companions. Consult those directly for usage details; they aren't duplicated here.

## Dependencies

Most workflows require **[CanlabCore](https://github.com/canlab/CanlabCore)** and **[SPM12](https://www.fil.ion.ucl.ac.uk/spm/)** on the MATLAB path. Neither is vendored into this repo — both must be cloned/installed separately and added to the path (see `prep/LaBGAScore_prep_s0_define_directories.m` and `pet/scripts/LaBGAScore_pet_a_set_up_paths_always_run_first.m`, which check for and set up these paths per study).

Several domain folders assume their own external toolbox, also not vendored, expected pre-installed and on the path:

| Folder | External toolbox |
|---|---|
| `cosmomvpa/` | [CoSMoMVPA](https://www.cosmomvpa.org/) |
| `decoding_toolbox/` | [The Decoding Toolbox (TDT)](https://sites.google.com/site/tdtdecodingtoolbox/) |
| `graphvar/` | [GraphVar](http://rfmri.org/GraphVar) |
| `juspace/` | [JuSpace](https://github.com/juryxy/JuSpace) |
| `mrs/` | [Osprey](https://github.com/schorschinho/osprey) |

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

**`firstlevel/`** — SPM/CANlab first-level GLM pipeline: `s1_options_dsgn_struct` (builds the CANlab-style `DSGN` design struct) → `s2_fit_model` (fits and diagnoses first-level models, cloning CANlab dependencies as needed) → `s3_diagnose_model` (publishes an HTML diagnostic report via CANlab's `scn_spm_design_check`, saves VIFs). A parallel phMRI-specific variant exists (`s1a`/`s1b`/`s2a`/`s2b`/`s3b`/`s3c`) for multisession/multitask and pharmacological-challenge designs. `functions/canlab_glm_subject_levels*_old.m` are older CANlab GLM functions kept because current example scripts still call them.

**`secondlevel/`** — group-level statistics and MVPA/ML pipelines: TFCE permutation inference (`group_tfce_from_subject_maps.m`, `tfce_one_fmri_dat.m`, `tfce_transform_3d.m`, `call_pTFCE.m`), PLS-DA/PLSR/Elastic Net pipelines with matching diagnostic-plotting functions (see the six `README_*.md` guides above), dice-overlap tools (`dice_statistic_image*.m`), an ROI/parcel extraction script (`LaBGAScore_secondlevel_extractparcels_sessions.m`), a MACS-toolbox model-space batch-setup script (`LaBGAScore_secondlevel_MS_mat_pipeline.m`), an MVPA-regression-on-connectivity-betas script (`LaBGAScore_secondlevel_mvpa_beta_maps_conn.m`), a PLS/ENet ROI-pipeline wrapper (`LaBGAScore_secondlevel_roi_run_plot_PLS_ENet_pipeline.m`), an object-oriented ML toolkit example (`LaBGAScore_secondlevel_ooFmriDataObjML_example.m`), and the `ProgressTracker` class.

**`stats_tools/`** — `LaBGAScore_Storey_FDR.m`, implementing Storey's positive FDR correction (falling back to Benjamini-Hochberg when its precondition isn't met).

**`atlas_mask_tools/`** — `LaBGAScore_atlas_binary_mask_from_atlas.m` and `LaBGAScore_atlas_rois_from_atlas.m` generate custom atlas/mask and per-ROI objects from a chosen atlas; the folder also ships a set of ready-made brain/gray-matter mask and template NIfTIs (`brain_masks/`, `brain_templates/`, `gray_matter_masks/`).

**`pet/`** — PET-specific pipeline: path setup (`a_set_up_paths_always_run_first.m`, `a2_set_default_options.m`), DICOM-to-BIDS conversion, preprocessing, kinetic modeling (e.g. `LaBGAScore_pet_model_TSPO_DPA714.m`), prep/signature-application scripts mirroring the fMRI `prep_*` naming, and a PLS/ENet pipeline runner. `functions/LCN_*.m` are the vendored legacy KU Leuven PET-processing functions.

**`mrs/`** — MRS pipeline built on Osprey: BIDS conversion (`LaBGAScore_prep_mrs2bids.m`), GE/Philips Osprey jobfiles (single- and multi-session), and run scripts (`LaBGAScore_mrs_run_osprey_GE.m`/`_Philips.m`).

**`cosmomvpa/`** — `LaBGAScore_cosmomvpa_searchlight_rsa.m`: first-level representational similarity analysis (behavioral vs. neural dissimilarity, leave-one-run-out cross-validation) plus group-level permutation/TFCE testing, via CoSMoMVPA.

**`decoding_toolbox/`** — TDT-based decoding: `LaBGAScore_decoding_template_xclass_acc.m` (first-level (cross-)classification accuracy plus group-level testing) and `LaBGAScore_decoding_SVM_between_subjects.m` (a full between-subject SVM decoding pipeline with permutation-based TFCE inference).

**`graphvar/`** — `LaBGAScore_prep_graphvar_input_from_conn.m`: builds GraphVar input files from CONN toolbox ROI-to-ROI connectivity output.

**`juspace/`** — `LaBGAScore_prep_juspace_input.m` (builds first-level con-image input for JuSpace PET-receptor spatial correlation analysis) and `LaBGAScore_juspace_corr_behav.m` (correlates behavioral data with JuSpace's spatial-correlation results).

**`power/`** — standalone power-analysis helper functions: `holmthreshold.m`, `medianIQR_to_meanSD.m`, `print_LEAR_matrix.m`, `sd_from_ci.m`.

**`figures/`** — plotting helpers: `canlabCmap.m` (CANlab colormap), `cluster_surface_plots.m`, `save_all_open_figures_smart.m`.

**`clean/`** — utilities: `LaBGAScore_clean_gzip_all_nii.m` (recursive `.nii` gzip), `LaBGAScore_clean_sourcedata.m` (cleans the `sourcedata` subdataset), `LaBGAScore_move_repos_matlabpath.m`, `LaBGAScore_smart_parallel_pool_setup.m`.

**`qr_code/`** — standalone Python utilities (`emailQR.py`, `email_QR_Input.py`, `QRtoPDF.py`) for generating and emailing QR codes; not MATLAB, and unrelated to the neuroimaging pipeline proper.

## Relationship to `CANlab_help_examples` (LaBGAS fork)

`CANlab_help_examples` (local path `/data/master_github_repos/CANlab_help_examples`, upstream [github.com/labgas/CANlab_help_examples](https://github.com/labgas/CANlab_help_examples)) is a sibling repo whose `Second_level_analysis_template_scripts/` directory provides LaBGAS's group-level (second-level) fMRI analysis templates — univariate GLM, cross-validated SVM, CANlab signature-pattern responses, searchlight correlation, and single-trial/runwise MVPA and mediation. Several of those scripts depend directly on LaBGAScore:

| LaBGAScore function/script | Called from (in CANlab_help_examples) | Purpose |
|---|---|---|
| `LaBGAScore_prep_s0_define_directories` (`prep/`) | `a_set_up_paths_always_run_first.m` | Defines `rootdir`, `codedir`, `BIDSdir`, `spmrootdir`, etc. for the study |
| `LaBGAScore_firstlevel_s1_options_dsgn_struct` (`firstlevel/`) | `a_set_up_paths_always_run_first.m` | Builds the first-level `DSGN` struct the second-level framework reads `DSGN.modeldir`/`DSGN.conditions`/`DSGN.contrastnames` from |
| `LaBGAScore_firstlevel_s2_fit_model.m` (`firstlevel/`) | upstream of `prep_3f_create_fmri_data_single_trial_object.m` (not called directly) | Fits first-level models; produces the single-trial con images that single-trial/runwise MVPA and mediation analyses consume |
| `LaBGAScore_atlas_binary_mask_from_atlas.m` (`atlas_mask_tools/`) | referenced via `atlasname_glm`/`atlasname_svm` options in `a2_set_default_options.m` | Generates custom `.mat` atlas/mask objects usable as GLM/SVM masks |
| `LaBGAScore_atlas_rois_from_atlas.m` (`atlas_mask_tools/`) | referenced via `roi_names`/`roi_modelname`/`roi_set_name` options | Generates per-ROI atlas objects for ROI-average analysis |

In short: LaBGAScore owns study setup, first-level modeling, and atlas/mask generation; `CANlab_help_examples` (LaBGAS fork) owns the second-level/group analysis templates built on top of LaBGAScore's outputs. For that repo's own internals, see its own [`README.md`](https://github.com/labgas/CANlab_help_examples/blob/main/Second_level_analysis_template_scripts/README.md)/`CLAUDE.md` under `Second_level_analysis_template_scripts/`.

## Tests and CI

There is no automated test suite and no CI pipeline in this repo. Verify changes manually by running the affected script against real or study data; several scripts support MATLAB's `publish()` to generate a date-stamped HTML report of their output (e.g. `cosmomvpa/LaBGAScore_cosmomvpa_searchlight_rsa.m`, `decoding_toolbox/LaBGAScore_decoding_template_xclass_acc.m`).

`clean/LaBGAScore_check_all_scripts.m` runs MATLAB's built-in Code Analyzer (`checkcode`) across every `.m` file in the repo (or a subtree you pass it) and reports the results, with genuine syntax errors called out separately from style/performance suggestions. It reliably catches parse errors — a file that cannot run past the flagged line — but does **not** catch undefined variables used at runtime, calls to functions that don't exist or aren't on the path, or logic bugs; those still require reading the code. Run it before committing as a fast baseline check, not a substitute for review.

## License

GNU General Public License v3.0 — see [`LICENSE`](LICENSE).
