# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

LaBGAScore holds LaBGAS's (Laboratory for Brain-Gut Axis Studies, KU Leuven) core/template MATLAB scripts for its standard neuroimaging analysis workflow. Scripts are either run directly from this repo, or (more commonly, for a real study) downloaded and adapted into that study's own `code` subdataset — every LaBGAS project is a DataLad superdataset (`proj_xxx`) with subdatasets `sourcedata`, `BIDS`, `derivatives`, `code`, `firstlevel`, `secondlevel`. GPLv3. No build system, test suite, or CI; this is a curated collection of runnable/copyable `.m` scripts, not a packaged product. One lightweight exception: `clean/LaBGAScore_check_all_scripts.m` (see "Static analysis" below).

Not everything here is MATLAB: `stats_tools/sas_macros/` holds SAS macros (`.sas`) for mixed-model effect sizes, and `qr_code/` holds standalone Python. Both are self-contained and outside the MATLAB conventions described below.

## Dependencies (not vendored)

Most workflows require **CANlabCore** and **SPM12** on the MATLAB path. Some domain folders assume their own external toolbox, also not vendored: `cosmomvpa/` → CoSMoMVPA, `decoding_toolbox/` → The Decoding Toolbox (TDT), `graphvar/` → GraphVar, `juspace/` → JuSpace, `mrs/` → Osprey, `stats_tools/sas_macros/` → SAS with SAS/STAT (no MATLAB involved). Exception: `pet/functions/LCN_*.m` are legacy KU Leuven PET-processing functions vendored directly into the repo. No package manager/manifest.

## Repository structure

Topic-organized top level, not a conventional toolbox layout: `prep/`, `firstlevel/`, `secondlevel/`, `stats_tools/`, `atlas_mask_tools/`, `pet/`, `mrs/`, `decoding_toolbox/`, `cosmomvpa/`, `graphvar/`, `juspace/`, `power/`, `figures/`, `clean/`, `qr_code/` (Python, not MATLAB). Most domains split into `<domain>/scripts/` + `<domain>/functions/`; `stats_tools/` splits into `functions/` (MATLAB) + `sas_macros/` (SAS). One class in the whole repo: `secondlevel/classes/ProgressTracker.m`. `secondlevel/` also has seven standalone usage guides (`README_ENet_*.md`, `README_PLSDA_*.md`, `README_PLSR_*.md`), `firstlevel/` has two (`README.md` for the task-fMRI chain, `README_phMRI.md` for the pharmacological-challenge chain), and `stats_tools/sas_macros/` and `clean/` have one each — all authoritative for their subject; don't duplicate their content.

## Script conventions

- Numbered/lettered prefixes (`prep_1_`, `s0_`/`s1_`/`s2_`/`s3_`, `a_`/`a2_`, and lettered variants `s1a_`/`s1b_`/`s2a_`/`s2b_`) mark ordered pipeline steps meant to be copied into a study's `code` subdataset and adapted — not generic library code.
- A `_example` suffix marks a script that is illustrative only.
- **Scripts** (no `function`/`classdef` declaration) use a comment header: `%% scriptname.m` title, `*USAGE*`, optionally `*OPTIONS*`/`*DEPENDENCIES*`/`*NOTES*`, then an author/date/version block — see `prep/LaBGAScore_prep_s0_define_directories.m`. **Functions** (`*/functions/` files, plus a few top-level ones like `power/holmthreshold.m`) use MATLAB's standard function help-text convention (H1 line, syntax) — a different, already-consistent convention, out of scope for the header-consistency work below.
- The `.sas` files follow neither: they use their own header block and the conventions recorded at the end of `stats_tools/sas_macros/README.md` (validate inputs and abort with `ERROR:` rather than emit a silently wrong dataset; state unverifiable assumptions in the log with `NOTE:`; suffix approximations `*_APPROX`).

## Relationship to sibling repo `CANlab_help_examples` (LaBGAS fork)

`CANlab_help_examples` (`/data/master_github_repos/CANlab_help_examples`, upstream `github.com/labgas/CANlab_help_examples`) provides LaBGAS's second-level (group) fMRI analysis templates, built on top of LaBGAScore:

| LaBGAScore function/script | Called from (in CANlab_help_examples) | Purpose |
|---|---|---|
| `LaBGAScore_prep_s0_define_directories` (`prep/`) | `a_set_up_paths_always_run_first.m` | Defines `rootdir`, `codedir`, `BIDSdir`, `spmrootdir`, etc. |
| `LaBGAScore_firstlevel_s1_options_dsgn_struct` (`firstlevel/`) | `a_set_up_paths_always_run_first.m` | Builds the first-level `DSGN` struct the second-level framework reads from |
| `LaBGAScore_firstlevel_s2_fit_model.m` (`firstlevel/`) | upstream of `prep_3f_create_fmri_data_single_trial_object.m` | Fits first-level models; produces single-trial con images |
| `LaBGAScore_atlas_binary_mask_from_atlas.m` (`atlas_mask_tools/`) | `atlasname_glm`/`atlasname_svm` options | Generates custom atlas/mask objects |
| `LaBGAScore_atlas_rois_from_atlas.m` (`atlas_mask_tools/`) | `roi_names`/`roi_modelname`/`roi_set_name` options | Generates per-ROI atlas objects |
| `LaBGAScore_smart_parallel_pool_setup.m` (`clean/`) | `c2a_second_level_regression.m`, `prep_3a_run_second_level_regression_and_save.m`, `prep_3c_run_SVMs_on_contrasts_masked.m` | Sets up the parallel pool before bootstrapping/permutation |
| `group_tfce_from_subject_maps.m` (`secondlevel/functions/`) | `prep_3a_run_second_level_regression_and_save.m` | Group TFCE from subject-level maps |
| `thresholded_fmri_data_from_statistic_image.m` (`secondlevel/functions/`) | `prep_3a_run_second_level_regression_and_save.m`, `c2_SVM_contrasts_masked.m` | Thresholded `fmri_data` object from a `statistic_image` |

The first two rows above are call-graph verified; the two `atlas_mask_tools` entries are
reached through option strings (`atlasname_glm`, `roi_names`) rather than direct calls, so
they do not appear in `dependencies.tsv`. The last three rows were found by the dependency
tooling and were previously undocumented. `DEPENDENCIES.md` in each repo is the generated,
authoritative version of this table.

LaBGAScore = study setup + first-level + atlas/mask generation; CANlab_help_examples (LaBGAS fork) = second-level templates built on top. See that repo's own `README.md`/`CLAUDE.md` under `Second_level_analysis_template_scripts/` for its internals.

## Documentation & audit history

Three documentation goals, mirroring work already done for `CANlab_help_examples` — all now complete:

1. **Rewrite `README.md`.** ✅ Done (commit `6edd5c5`). Replaced the 3-line stub with an extensive document covering: what this repo is; the DataLad project structure; repository structure/naming conventions; dependencies; a domain-by-domain overview (one paragraph per top-level folder, its entry-point script(s), pointers to richer docs like `secondlevel/README_*.md`); the `CANlab_help_examples` dependency table above; the no-tests/CI note; license. Modeled on `CANlab_help_examples/Second_level_analysis_template_scripts/README.md`'s structure (TOC with anchors).

2. **Normalize script header comments — structural pass.** ✅ Done (commit `654f450`). Scope: the 47 true MATLAB **script** files repo-wide (no `function`/`classdef` declaration — see "Script conventions" above; function files and `ProgressTracker.m` are out of scope). Standardized separator style (dashes, not the old `%____...` underscore block), section-label formatting (asterisked `*USAGE*`/`*OPTIONS*`/`*DEPENDENCIES*`/`*NOTES*`, not plain), and dropped the old `@(#)%` SCCS-style prefix on the version-stamp line. This was a structural/formatting pass only — it preserved existing substantive content and did not verify that content was accurate or complete.

3. **Review script header content against actual code — accuracy pass.** ✅ Done (commits `27d8d8a`, `f53dacb`). For the same 47 script files, read each script's full body and checked whether its structurally-normalized header content was accurate, complete, and current. This also surfaced 19 real code defects (undefined variables, wrong/nonexistent function calls, hardcoded study-specific setup calls left in generic templates, copy-paste errors) unrelated to documentation, all fixed alongside the header content.

## Provenance and dependency tooling

`clean/LaBGAScore_prov_*` and `clean/LaBGAScore_dep_*` record which version of each
dependency produced a given result, and generate the dependency documentation. Full
detail in `clean/README_provenance.md`; the essentials:

- **The dependencies are the reproducibility gap.** Scripts are frozen per project in the
  `code` subdataset; CanlabCore and the other clones under `/data/master_github_repos` are
  not. Vendoring them per project was considered and rejected (Neuroimaging_Pattern_Masks
  is 18 GB, and pinning blocks bug fixes).
- **`LaBGAScore_prov_publish`** replaces the bare `publish(...)` call documented in each
  script header, and adds a provenance table to the html report plus a `.tsv`/`.mat` in
  `results/notes/`. No script templates were edited — deliberately, since every study has
  its own renamed copies.
- **`LaBGAScore_prov_resolve_retrospective`** reconstructs the same information for
  analyses already run, from each artifact's own embedded date and each clone's git
  reflog. Writes sidecars; never modifies existing files.
  - It covers **result `.mat` files as well as published reports** — only a few scripts
    per model are ever published, so reports alone leave most of the pipeline
    undocumented (in `proj_cfs`: 29 reports against 98 result files). A `.mat` carries a
    `Created on:` stamp in its 116-byte header, which dates it the way a report's
    `DC.date` does, with a time of day and surviving git-annex.
  - Output is **one page per run**, not per file: a script's report and its `.mat` files
    are grouped together, ranked by evidence quality. The group key includes the subject
    where there is one, since at first level the same script runs once per subject.
  - It follows whichever directory layout a model already has — `results/notes` and
    `results/html` at second level, `<model>/provenance/` at first level, where no
    `results/` directory exists.
- **Reflogs are the critical, perishable evidence.** `gc.reflogExpire` defaults to 90 days
  and prunes silently. `clean/labgascore_prov_protect_reflogs.sh` disables that; it has
  been run on all 29 clones on the LaBGAS server.
- **Nothing shells out.** Commit, branch, remote, reflog and dirty state are read straight
  from the plain-text files under `.git`. `LaBGAScore_prov_gitstatus.m` reimplements git's
  own dirty check and is verified to reproduce `git status --porcelain` exactly.
- **`DEPENDENCIES.md`, `dependencies.tsv` and `dependencies.yml` are generated** by
  `clean/LaBGAScore_dep_report.m`. Never hand-edit them. The LaBGAS website collects the
  `.yml` via `scripts/refresh_dependencies.py`. In this repo they sit at the root and cover
  all ~130 files; in `CANlab_help_examples` they sit in
  `Second_level_analysis_template_scripts/` and cover exactly the 19 scripts listed in that
  folder's README, not the ~113 scripts present.
- **`LaBGAScore_check_display`** reports whether the current X2go session can produce a
  full-size report figure, and what to change if not. `publish()` captures figures from the
  screen, so session size and DPI decide how figures come out; recommended settings per
  screen are in `LaBGAS_fMRI_analysis_workflow.md` ("Set up your X2go display for
  publishing figures"). `LaBGAScore_prov_publish` records the screen geometry, DPI and the
  resulting figure dimensions in every report, and flags figures whose size was set by the
  display rather than by the script.

Resolution is index-based rather than `which()`-based, because `which` returns one
path-order-dependent answer and on this setup it is the wrong one for the calls that
matter most: `which('predict')` returns a liblinear `.mexa64` from SPM's decoding
toolbox and `which('ttest')` returns MATLAB's Statistics Toolbox, where the second-level
scripts actually call `@fmri_data/predict` and `@fmri_data/ttest`. `which()` is used only
to break ties among candidates the index already found.

Six rules stop the call graph inventing dependencies (classdef property declarations are
not calls; dot-calls and ambiguous names may reinforce a repository but never introduce
one; core-language builtins are never third-party calls; evidence must be repo-unique;
genuine cross-repo ties are broken with `which()`). Together they cut the `proj_cfs`
second-level record from 1639 rows to 552 and removed BrainSpace, gift, cocoanCORE,
ExploreASL and others that nothing here calls. See `clean/README_provenance.md` — do not
relax these without re-checking the negative controls documented there.

**State as of 2026-09-01:** the retrospective has been run over `proj_cfs` and
`proj_discoverie` (second and first level). Those outputs are written but **not committed**
in the `proj_*` datasets, pending review.

## Static analysis

`clean/LaBGAScore_check_all_scripts.m` runs MATLAB's built-in Code Analyzer (`checkcode`) across every `.m` file in the repo (or a subtree passed as its `rootdir` argument) and prints a report, separating genuine syntax errors from style/performance suggestions. Calibrated against real bugs found during the accuracy pass above: `checkcode` reliably catches parse errors (e.g. it did catch the `sort{}` syntax bug fixed in commit `27d8d8a`) but does **not** catch undefined variables used at runtime, calls to functions that don't exist or aren't on the path, or logic bugs (e.g. wrong array indexing) — those all had to be found by reading the code, not by static analysis. Treat it as a fast baseline check run before committing, not a substitute for the kind of full read-through that found the 19 defects above.

`matlab.codetools.requiredFilesAndProducts` is **not** usable in this tree and should not
be re-attempted: it aborts entirely on a syntax error anywhere in the transitive closure,
and there are 33 such files (CanlabCore 10, CanlabPrivate 21, canlab_single_trials 1,
CANlab_help_examples 1). `clean/LaBGAScore_dep_map.m` records those as `unparseable` and
carries on.

It walks `.m` files only, so `stats_tools/sas_macros/` is covered by nothing. Those macros have also not yet been executed against a real model in SAS — their arithmetic was verified independently and their syntax checked by reading, no more. That caveat is stated at the top of their README and should be removed there once someone has run them.
