# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

LaBGAScore holds LaBGAS's (Laboratory for Brain-Gut Axis Studies, KU Leuven) core/template MATLAB scripts for its standard neuroimaging analysis workflow. Scripts are either run directly from this repo, or (more commonly, for a real study) downloaded and adapted into that study's own `code` subdataset — every LaBGAS project is a DataLad superdataset (`proj_xxx`) with subdatasets `sourcedata`, `BIDS`, `derivatives`, `code`, `firstlevel`, `secondlevel`. GPLv3. No build system, linter, test suite, or CI — none should be added; this is a curated collection of runnable/copyable `.m` scripts, not a packaged product.

## Dependencies (not vendored)

Most workflows require **CANlabCore** and **SPM12** on the MATLAB path. Some domain folders assume their own external toolbox, also not vendored: `cosmomvpa/` → CoSMoMVPA, `decoding_toolbox/` → The Decoding Toolbox (TDT), `graphvar/` → GraphVar, `juspace/` → JuSpace, `mrs/` → Osprey. Exception: `pet/functions/LCN_*.m` are legacy KU Leuven PET-processing functions vendored directly into the repo. No package manager/manifest.

## Repository structure

Topic-organized top level, not a conventional toolbox layout: `prep/`, `firstlevel/`, `secondlevel/`, `stats_tools/`, `atlas_mask_tools/`, `pet/`, `mrs/`, `decoding_toolbox/`, `cosmomvpa/`, `graphvar/`, `juspace/`, `power/`, `figures/`, `clean/`, `qr_code/` (Python, not MATLAB). Most domains split into `<domain>/scripts/` + `<domain>/functions/`. One class in the whole repo: `secondlevel/classes/ProgressTracker.m`. `secondlevel/` also has six standalone usage guides (`README_ENet_*.md`, `README_PLSDA_*.md`, `README_PLSR_*.md`) — don't duplicate their content.

## Script conventions

- Numbered/lettered prefixes (`prep_1_`, `s0_`/`s1_`/`s2_`/`s3_`, `a_`/`a2_`, and lettered variants `s1a_`/`s1b_`/`s2a_`/`s2b_`) mark ordered pipeline steps meant to be copied into a study's `code` subdataset and adapted — not generic library code.
- A `_example` suffix marks a script that is illustrative only.
- **Scripts** (no `function`/`classdef` declaration) use a comment header: `%% scriptname.m` title, `*USAGE*`, optionally `*OPTIONS*`/`*DEPENDENCIES*`/`*NOTES*`, then an author/date/version block — see `prep/LaBGAScore_prep_s0_define_directories.m`. **Functions** (`*/functions/` files, plus a few top-level ones like `power/holmthreshold.m`) use MATLAB's standard function help-text convention (H1 line, syntax) — a different, already-consistent convention, out of scope for the header-consistency work below.

## Relationship to sibling repo `CANlab_help_examples` (LaBGAS fork)

`CANlab_help_examples` (`/data/master_github_repos/CANlab_help_examples`, upstream `github.com/labgas/CANlab_help_examples`) provides LaBGAS's second-level (group) fMRI analysis templates, built on top of LaBGAScore:

| LaBGAScore function/script | Called from (in CANlab_help_examples) | Purpose |
|---|---|---|
| `LaBGAScore_prep_s0_define_directories` (`prep/`) | `a_set_up_paths_always_run_first.m` | Defines `rootdir`, `codedir`, `BIDSdir`, `spmrootdir`, etc. |
| `LaBGAScore_firstlevel_s1_options_dsgn_struct` (`firstlevel/`) | `a_set_up_paths_always_run_first.m` | Builds the first-level `DSGN` struct the second-level framework reads from |
| `LaBGAScore_firstlevel_s2_fit_model.m` (`firstlevel/`) | upstream of `prep_3f_create_fmri_data_single_trial_object.m` | Fits first-level models; produces single-trial con images |
| `LaBGAScore_atlas_binary_mask_from_atlas.m` (`atlas_mask_tools/`) | `atlasname_glm`/`atlasname_svm` options | Generates custom atlas/mask objects |
| `LaBGAScore_atlas_rois_from_atlas.m` (`atlas_mask_tools/`) | `roi_names`/`roi_modelname`/`roi_set_name` options | Generates per-ROI atlas objects |

LaBGAScore = study setup + first-level + atlas/mask generation; CANlab_help_examples (LaBGAS fork) = second-level templates built on top. See that repo's own `README.md`/`CLAUDE.md` under `Second_level_analysis_template_scripts/` for its internals.

## Current objectives in this repo

Two active documentation goals, mirroring work already done for `CANlab_help_examples`:

1. **Rewrite `README.md`.** It's currently a 3-line stub. Replace it with an extensive document covering: what this repo is; the DataLad project structure; repository structure/naming conventions; dependencies; a domain-by-domain overview (one paragraph per top-level folder, its entry-point script(s), pointers to richer docs like `secondlevel/README_*.md`); the `CANlab_help_examples` dependency table above; the no-tests/CI note; license. Model it on `CANlab_help_examples/Second_level_analysis_template_scripts/README.md`'s structure (TOC with anchors).

2. **Normalize script header comments.** Scope: the 47 true MATLAB **script** files repo-wide (no `function`/`classdef` declaration — see "Script conventions" above; function files and `ProgressTracker.m` are out of scope). Most already approximate the target header shape but are inconsistent in details: separator style (old `%____...` underscore block vs. `% ----...` dashes — standardize on dashes), section-label formatting (plain `USAGE` vs. asterisked `*USAGE*` — standardize on asterisked), author/date block wording, and an old `@(#)%` SCCS-style prefix on the version-stamp line (drop it). A few scripts (`decoding_toolbox/LaBGAScore_decoding_SVM_between_subjects.m`, `clean/LaBGAScore_smart_parallel_pool_setup.m`) lack the structured header entirely and need one built from scratch. This is a structural/formatting pass — preserve existing substantive content (authors, dates, version numbers, USAGE/NOTES prose); only add new content (e.g. an `*OPTIONS*` section) where a script has clearly user-configurable variables not yet documented, derived from reading the script body, not invented.
