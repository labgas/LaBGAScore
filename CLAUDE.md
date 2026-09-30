# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

LaBGAScore holds LaBGAS's (Laboratory for Brain-Gut Axis Studies, KU Leuven) core/template MATLAB scripts for its standard neuroimaging analysis workflow. Scripts are either run directly from this repo, or (more commonly, for a real study) downloaded and adapted into that study's own `code` subdataset — every LaBGAS project is a DataLad superdataset (`proj_xxx`) with subdatasets `sourcedata`, `BIDS`, `derivatives`, `code`, `firstlevel`, `secondlevel`. GPLv3. No build system, test suite, or CI; this is a curated collection of runnable/copyable `.m` scripts, not a packaged product. One lightweight exception: `clean/LaBGAScore_check_all_scripts.m` (see "Static analysis" below).

Not everything here is MATLAB: `stats_tools/sas_macros/` holds SAS macros (`.sas`) for mixed-model effect sizes, and `qr_code/` holds standalone Python. Both are self-contained and outside the MATLAB conventions described below.

## Dependencies (not vendored)

Most workflows require **CANlabCore** and **SPM12** on the MATLAB path. Some domain folders assume their own external toolbox, also not vendored: `cosmomvpa/` → CoSMoMVPA, `decoding_toolbox/` → The Decoding Toolbox (TDT), `graphvar/` → GraphVar, `juspace/` → JuSpace, `mrs/` → Osprey, `stats_tools/sas_macros/` → SAS with SAS/STAT (no MATLAB involved). Exception: `pet/functions/LCN_*.m` are legacy KU Leuven PET-processing functions vendored directly into the repo. No package manager/manifest.

## Repository structure

Topic-organized top level, not a conventional toolbox layout: `prep/`, `firstlevel/`, `secondlevel/`, `stats_tools/`, `atlas_mask_tools/`, `pet/`, `mrs/`, `decoding_toolbox/`, `cosmomvpa/`, `graphvar/`, `juspace/`, `power/`, `figures/`, `clean/`, `qr_code/` (Python, not MATLAB). Most domains split into `<domain>/scripts/` + `<domain>/functions/`; `stats_tools/` splits into `functions/` (MATLAB) + `sas_macros/` (SAS). One class in the whole repo: `secondlevel/classes/ProgressTracker.m`. `secondlevel/` has a `README.md` indexing the folder plus seven standalone usage guides (`README_ENet_*.md`, `README_PLSDA_*.md`, `README_PLSR_*.md`), `firstlevel/` has two (`README.md` for the task-fMRI chain, `README_phMRI.md` for the pharmacological-challenge chain), and `prep/`, `decoding_toolbox/`, `atlas_mask_tools/` and `stats_tools/sas_macros/` have one each, while `clean/` has two (`README.md` indexing the folder, `README_provenance.md` authoritative for the provenance and dependency tooling) — all authoritative for their subject; don't duplicate their content.

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
| `tfce_fwe_from_null.m` (`secondlevel/functions/`) | `c2a_second_level_regression.m` | Max-statistic FWE p-values from a saved TFCE permutation null |
| `LaBGAScore_region_table.m` (`secondlevel/functions/`) | `c2a_second_level_regression.m` | `@region/table` vendored without the large-value clipping. Upstream ends `get_signed_max` with `maxZ = norminv(1 - 1E-12)` (= **7.0345**) and clips every value above it, not only the infinities its own comment describes. As of c2a v8.5 **all five** table branches use this copy; before that only the two TFCE branches did, and FDR / uncorrected / Bayesian peaks above 7.0345 were silently flattened. **Also relabels the Bayes column (2026-09-25).** `@statistic_image/estimateBayesFactor` ends with `BF.dat = 2*log(bf10)` and sets `.type='BF'`, so the column was headed `maxBF` while holding **2·ln(BF10)** - a printed 15.05 is a Bayes factor of 1854, not 15. The column is now `max_2lnBF` with `maxBF10 = exp(max_2lnBF/2)` beside it. Note how the two defects compounded: clipping at 7.0345 on a 2·ln scale caps the reported evidence at BF10 ≈ 34, so an unclipped peak of 39591 printed as 7.03 |
| `LaBGAScore_region_table_safe.m` (`secondlevel/functions/`) | `c2a_second_level_regression.m` | Three-output wrapper so an undisplayable contrast does not abort the report |

The first two rows above are call-graph verified; the two `atlas_mask_tools` entries are
reached through option strings (`atlasname_glm`, `roi_names`) rather than direct calls, so
they do not appear in `dependencies.tsv`. The last three rows were found by the dependency
tooling and were previously undocumented. `DEPENDENCIES.md` in each repo is the generated,
authoritative version of this table.

LaBGAScore = study setup + first-level + atlas/mask generation; CANlab_help_examples (LaBGAS fork) = second-level templates built on top. See that repo's own `README.md`/`CLAUDE.md` under `Second_level_analysis_template_scripts/` for its internals.

## Bayes factors: the stored scale is 2·ln(BF)

`@statistic_image/estimateBayesFactor` returns **2·ln(BF10)** (Kass & Raftery 1995),
positive favouring H1. Three consequences worth not rediscovering:

- **Converting.** `BF10 = exp(value/2)`, never `exp(value)`. Getting this wrong
  inflates the apparent evidence *for the null* dramatically: on model_2c_IOM it
  turned "83% moderate evidence for the null, median BF 0.16" into a spurious
  "82% strong evidence for the null, median BF 0.027".
- **Thresholds in the pipeline are correct.** `c2a` converts on use with
  `2*log(BF_threshold_glm)`; `prep_3a` hardcodes `2.1972`, which is `2*ln(3)` and
  is labelled "|BF| > 3". The bare constant is what misleads - it is not ln(9).
- **The JZS BF has a sample-size floor.** `t1smpbf(0, n)` bounds how much evidence
  for H0 is attainable: n=70 floors at BF10 = 0.131 (7.6:1), n=93 at 0.115,
  n=158 at 0.089. **Below roughly n=100, "strong evidence for the null"
  (BF10 < 1/10) is unreachable no matter how null the data are.** A Bayes map
  showing 0% strong-for-null in a small sample is reporting a design ceiling,
  not a weak result, and should say so.

## Documentation & audit history

Four documentation goals, mirroring work already done for `CANlab_help_examples` — all now complete:

1. **Rewrite `README.md`.** ✅ Done (commit `6edd5c5`). Replaced the 3-line stub with an extensive document covering: what this repo is; the DataLad project structure; repository structure/naming conventions; dependencies; a domain-by-domain overview (one paragraph per top-level folder, its entry-point script(s), pointers to richer docs like `secondlevel/README_*.md`); the `CANlab_help_examples` dependency table above; the no-tests/CI note; license. Modeled on `CANlab_help_examples/Second_level_analysis_template_scripts/README.md`'s structure (TOC with anchors).

2. **Normalize script header comments — structural pass.** ✅ Done (commit `654f450`). Scope: the 47 true MATLAB **script** files repo-wide (no `function`/`classdef` declaration — see "Script conventions" above; function files and `ProgressTracker.m` are out of scope). Standardized separator style (dashes, not the old `%____...` underscore block), section-label formatting (asterisked `*USAGE*`/`*OPTIONS*`/`*DEPENDENCIES*`/`*NOTES*`, not plain), and dropped the old `@(#)%` SCCS-style prefix on the version-stamp line. This was a structural/formatting pass only — it preserved existing substantive content and did not verify that content was accurate or complete.

3. **Review script header content against actual code — accuracy pass.** ✅ Done (commits `27d8d8a`, `f53dacb`). For the same 47 script files, read each script's full body and checked whether its structurally-normalized header content was accurate, complete, and current. This also surfaced 19 real code defects (undefined variables, wrong/nonexistent function calls, hardcoded study-specific setup calls left in generic templates, copy-paste errors) unrelated to documentation, all fixed alongside the header content.

4. **Per-domain deep pass.** ✅ Done for `prep/` (commit `0976ed2`), `decoding_toolbox/` (commits `1e124b8`, `dc9e313`), `atlas_mask_tools/`, `clean/` and `secondlevel/`. Goals 2 and 3 were repo-wide but per-file; this one takes a single domain and documents it as a whole. Added `prep/README.md` (workflow-position diagram, generic-vs-example status per script, the variables `prep_s0` puts in the base workspace and who consumes them, the `events.tsv` contract with `firstlevel`, known limitations) and expanded the six `prep/` headers to the `firstlevel` style — numbered step-by-step `*USAGE*` with `NOTE:` lines naming the functions each step calls. Reading the bodies for that surfaced **six more real defects** beyond the 19 from goal 3, all fixed: `pheno_tsv = false` crashed both `s1` scripts on an undefined `pheno_file`; the duplicate-logfile check tested `size(logfilepath,1) > 1` on a single path string and so could never fire, silently skipping the run as "logfile missing"; single-session `s1` referenced an undefined `subjs{sub}`; an empty `time_zero` was uncaught and yielded an empty table; `gzip('s6*')` hardcoded the default prefix in both `s2` scripts, so changing `prefix` left the smoothed images unzipped for `delete('*.nii')` to remove — **silently destroying the smoothing output**; and `s0` derived `spmrootdir` by splitting on `'/spm.m'`, which fails on Windows. The two `s2` scripts were also deduplicated (the smoothing body existed in four verbatim copies) and made consistent with each other. **The lesson for the remaining domains:** goal 3's per-file pass did not catch any of these, because they are cross-branch and cross-script inconsistencies — an option guarded in one place but not another, two files that should agree and don't — which only a whole-domain read surfaces. Note also that `checkcode` reports 0 parse-level messages across `prep/` both before and after, and that these fixes are verified by reading only, not by running.

   `decoding_toolbox/` was the second domain, and it repeated the pattern with one new twist. Added `decoding_toolbox/README.md` (how the two pipelines differ, the four levers for a site/nuisance confound and the label-blindness invariant across them, which feature space each route decodes, per-script audit status) and rewrote the header of `LaBGAScore_decoding_template_xclass_acc.m`. The twist: that script **did not parse**. A `signrank` call was missing a closing parenthesis, which `checkcode` reports as `NOPAR`, so it could not have run in any form since its last use in 2023 — and goals 2 and 3 both read this file without noticing, because they were reading headers against code rather than running the checker. **Run `clean/LaBGAScore_check_all_scripts.m` on a domain before reading it.** The other fixes were of the same cross-cutting kind as `prep/`: `labelnames` indexed by condition number rather than position (so any non-contiguous `conds2include` silently produced empty labels), a subject count taken from a leaked loop counter, a mask resampled to subject 1 and reused unchecked, and a chance level that assumed TDT's confusion matrix is in percent — which, on a proportion-scaled matrix, would have returned a "significant" result for every condition. Several settled questions were also borrowed from the domain's own tested script (`scaling_regime`, `plot_selected_voxels`, `plugin_set_figure_size`) rather than re-derived, which is the cheapest part of a per-domain pass and only possible because both scripts were in view at once.

   `atlas_mask_tools/` was the third, and needed no code changes - the two scripts had already been audited and their headers were accurate. What it needed was the **mask library documented by measurement**. The naming conventions were undocumented and not guessable: for the gray-matter masks the trailing number is the GM probability ×100 applied to the bundled probseg (verified - the voxel count matches thresholding that map exactly at all six thresholds), while for the brain masks it is a **T1w intensity cut in the template's own 0-10000 units, applied within the template brain mask** - `_0` is byte-for-byte templateflow's `desc-brain_mask`, `_1000` is exactly that ∩ (T1w ≥ 1000), and `_1800` matches ∩ (T1w ≥ 1800) to within 2 voxels of 1.88 M. In both families higher means sparser, but the GM masks span 1557→1205 ml across the shipped range while the brain masks differ by under a quarter of one percent, which is worth knowing before anyone agonises over the choice. **The lesson: for a domain that ships data rather than only code, the documentation task is measurement, not reading.** Voxel counts, grids and nesting relations are facts about the files that no amount of header review would have produced.

   `clean/` and `secondlevel/` needed indexes rather than audits - both are recent work, and `secondlevel/` already had its own `CLAUDE.md` plus seven ML guides. Each new README says at the top what it does NOT cover, so authority stays in one place per subject. Writing the `secondlevel/` one nevertheless turned up two documentation defects: `secondlevel/CLAUDE.md` still headed its audit section *"findings, not yet fixed"* when every item in it is marked FIXED or RESOLVED (corrected), and the cross-repo dependency table in the root `README.md` was missing six of the eleven rows this file already carried - every `secondlevel/functions/` entry plus `smart_parallel_pool_setup` (added). The root README also framed LaBGAScore as owning "study setup, first-level modeling, and atlas/mask generation" while the fork owned the second level, which undersells `secondlevel/`: the TFCE stack, the ML pipelines and the reporting functions live HERE and the fork owns the templates that drive them. **Two tables covering the same relationship drift apart; the generated `DEPENDENCIES.md` is the one to trust.**

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
  `Second_level_analysis_template_scripts/` and cover exactly the 20 scripts listed in that
  folder's README, not the ~113 scripts present.
- **`LaBGAScore_check_display`** reports whether the current X2go session can produce a
  full-size report figure, and what to change if not. It applies to the **interactive**
  route only: there, `publish()` captures figures from the screen, so session size and DPI
  decide how figures come out. Recommended settings per screen are in
  `LaBGAS_fMRI_analysis_workflow.md` ("Setting up X2go, if you want higher-resolution
  figures"). `LaBGAScore_prov_publish` records the screen geometry, DPI and the resulting
  figure dimensions in every report, and flags figures whose size was set by the display
  rather than by the script.

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

**State as of 2026-09-23:** the retrospective has been run over `proj_cfs` and
`proj_discoverie`, second and first level. Commit status differs per dataset and was
re-checked on this date:

| dataset | provenance sidecars tracked | untracked |
|---|---|---|
| `proj_discoverie/secondlevel` | **46** | 0 |
| `proj_discoverie/firstlevel` | 0 | 332 |
| `proj_cfs/secondlevel` | 0 | 57 |
| `proj_cfs/firstlevel` | 0 | 137 |

So `proj_discoverie/secondlevel` is committed; the other three are still written-but-not-committed,
pending review. Re-derive these counts rather than trusting them — `git ls-files | grep -c provenance`
against `git status --porcelain | grep -c provenance` in each subdataset — since they move whenever
someone runs a `datalad save`.

## Running scripts and publishing reports

**Headless is the default**; X2go is the alternative when higher-resolution figures are
wanted. See "Running scripts and publishing reports" in `README.md` and section 2 of
`LaBGAS_fMRI_analysis_workflow.md`. Facts worth not rediscovering:

- **`matlab -batch` cannot `publish()`** — "Unable to run the 'publish' function, because
  it is not supported for this ...". `-nodisplay` with `-r` works. Verified directly.
- **Redirect stdin from `/dev/null`** with `-r`, or MATLAB can exit before running
  anything, which looks like a silent no-op.
- Headless reports `ScreenSize [1 1 1024 768]`, `ScreenPixelsPerInch 72`,
  `feature('ShowFigureWindows') == 0`. Figures are **not** capped by that virtual screen:
  without a display `publish` prints rather than screen-captures, so density grids come out
  at e.g. 1440x2160 px. The only cost of headless is pixel density (72 dpi against 96-144).
- **`publish()` catches a script error into the html and returns normally**, so a crashed
  run raises nothing and exits 0. This is the single biggest source of silently wrong
  results in this pipeline. `clean/LaBGAScore_run_reports.m` exists to close it: it reads
  each report back and fails on `<pre class="codeoutput error">`, and optionally asserts the
  script's output file exists. Match that **markup**, never the report text — phrases like
  "Error in" and "Unrecognized function" occur in ordinary comments, which `publish` renders
  as prose, so a text search flags scripts that ran perfectly. Verified both ways.
- MATLAB's `regexp` does **not** support `\b` as a word boundary (it uses `\<`/`\>`). A
  pattern like `\berror\b` silently never matches; this made the detector above report
  every failing report as OK until it was found by testing against a deliberately failing
  script.
- `clean/labgascore_run_headless.sh` wraps all of the above. Its MATLAB code goes into a
  temporary `.m` file rather than a one-line `-r` string, because flattening newlines makes
  a MATLAB `%` comment swallow the rest of the command and runs `for`/`if` headers into
  their bodies.

## Static analysis and the script checkers

`clean/LaBGAScore_check_all_scripts.m` runs MATLAB's built-in Code Analyzer (`checkcode`) across every `.m` file in the repo (or a subtree passed as its `rootdir` argument) and prints a report, separating genuine syntax errors from style/performance suggestions. Calibrated against real bugs found during the accuracy pass above: `checkcode` reliably catches parse errors (e.g. it did catch the `sort{}` syntax bug fixed in commit `27d8d8a`) but does **not** catch undefined variables used at runtime, calls to functions that don't exist or aren't on the path, or logic bugs (e.g. wrong array indexing) — those all had to be found by reading the code, not by static analysis. Treat it as a fast baseline check run before committing, not a substitute for the kind of full read-through that found the 19 defects above.

`checkcode` is blind to two *ordering* failures that cost hours each when they
happened, so `clean/` carries a Python checker for each. Run both on a study's
model script directory before launching any long chain.

- **`clean/use_before_def.py`** — an option **read above the line that defines
  it**. A guarded default only protects an option if it runs before every read;
  guard it beside one use and miss an earlier one and the script dies on
  *"Unrecognized function or variable"* at the earlier line. Real case:
  `contrast_objects_tag` guarded beside `savefilenamedata` while a `printhdr`
  eight lines above also read it — `prep_3` died there after 75 minutes, having
  saved nothing.
- **`clean/set_after_use.py`** — an option **set below the line that already
  consumed it**. Nothing errors; the default is simply what runs. Real case: a
  study copy set `results_tag` ~290 lines below the block that builds
  `tdt_resultsdir` from it, so eleven decoding runs wrote into the same untagged
  directory and overwrote each other's maps. Statistics were unaffected
  (`cfg.results.overwrite = 1`), the artefacts were not. **Advisory, not a
  gate** — a variable legitimately reused for successive outputs
  (`savefilename`, `figtitle`, `varnames`) trips it; expect a handful per model
  and read them rather than chase them.

Both ship positive controls in `clean/checker_positive_controls/`, run by
`run_controls.sh`, which asserts that each checker flags its own control and
ignores the other's. **`checkcode` reports zero messages on either control
file** — measured, not assumed — which is precisely why these exist alongside
it. If a change to a checker makes its control pass, the checker is broken.

These two automate traps 6 and 7 of the ten catalogued in
`LaBGAS_fMRI_analysis_workflow.md`; the rest still need reading.

`matlab.codetools.requiredFilesAndProducts` is **not** usable in this tree and should not
be re-attempted: it aborts entirely on a syntax error anywhere in the transitive closure,
and there are 33 such files (CanlabCore 10, CanlabPrivate 21, canlab_single_trials 1,
CANlab_help_examples 1). `clean/LaBGAScore_dep_map.m` records those as `unparseable` and
carries on.

It walks `.m` files only, so `stats_tools/sas_macros/` is covered by nothing. Those macros have also not yet been executed against a real model in SAS — their arithmetic was verified independently and their syntax checked by reading, no more. That caveat is stated at the top of their README and should be removed there once someone has run them.
