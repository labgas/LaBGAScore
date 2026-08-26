# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Scope

This file is scoped to `secondlevel/`. The repo-root `CLAUDE.md` (one level up) covers repo-wide
concerns — DataLad project layout, cross-repo relationship to `CANlab_help_examples`, naming
conventions, header-normalization audit history — and is not repeated here.

The seven `README_*.md` files in this folder are the authoritative user guides for the PLS-DA /
PLSR / Elastic Net pipelines and their plotting companions (options, inputs, every `results` field). Read
them before changing those functions; do not duplicate their content elsewhere.

### Completed: ML pipeline family overhaul

The four ML pipelines were reworked along three axes. The code and the seven `README_*.md` guides
are in sync as of this work.

1. **Fold-wise covariate residualization.** The pipelines previously had no way to control for
   covariates inside the CV loop — both calling scripts told the user to residualize beforehand,
   which fits the nuisance model on the full sample and so uses test-fold information in the
   nuisance betas. Replaced by train-only nuisance regression per fold (`opts.covariates`), at all five
   preprocessing sites. Permutation is Freedman-Lane on Y for PLSR (continuous outcome) and plain
   label permutation for the three classifiers, which is valid there because X is residualized
   in-fold and so is already orthogonal to the covariates (the Kennedy scheme). The paired pipeline
   keeps its within-subject label swap, which preserves subject-level covariates exactly.

   **The direction of the resulting bias is not what one might assume, and was measured rather than
   assumed.** Full-sample residualization does *not* simply inflate performance. When the covariate
   explains most of the apparent X–Y association, over-removal drives the estimate far *below*
   chance (measured: PLS-DA AUC 0.33 and PLSR Q² −0.77 on data where the true answer is chance,
   against 0.52 / −0.08 for the fold-wise route). When X carries genuine signal beyond the confound,
   the two routes agree closely (0.720 vs 0.731). This matches the known result that confound
   regression on a full dataset can produce systematically below-chance decoding. The case for the
   fold-wise version is therefore that it is *unbiased and leakage-free*, giving the right answer in
   both regimes — not that it removes optimism.
2. **De-duplication of the local subfunctions.** Done — see the helper table below. This was a
   prerequisite for (1), not a parallel cleanup: covariate handling has to reach five preprocessing
   sites per pipeline, and duplicated helpers would have meant twenty chances to miss one. The
   extraction was verified behaviour-preserving by diffing all 148 `results` fields across the four
   pipelines before and after, under a serial pool.
3. **Defect sweep.** Fixed, in rough order of severity:
   - *Permutation p-values were badly miscalibrated.* The null came from the untuned `quickCV_*`
     estimator while the observed statistic came from tuned repeated nested CV. Measured on
     true-null data this gave a **false-positive rate near 40 % at alpha = 0.05**. Both sides now use
     the same tuned estimator (`results.quickCV_observed`), with the `(sum+1)/(n+1)` convention.
     **p-values from earlier runs are inflated.**
   - *ENet hyperparameter selection* maximized a single inner fold's AUC over a flattened
     alpha × fold × lambda grid. It does not inflate the outer AUC (the outer fold is held out from
     selection), but it selected roughly 20x too little regularization — median lambda 0.068 vs 1.47
     — giving unstable, dense weight maps. Now fold-averaged with a single refit, via
     `selectENetHyperparams`.
   - *Inner folds inherited the outer fold's preprocessing constants* in three of four pipelines.
     Each inner fold now preprocesses itself.
   - *RNG was not reproducible under `parfor`.* Now seeded per iteration via `setParforStream`;
     verified identical across 2- and 4-worker pools. ENet's headline AUC was previously not
     reproducible at all.
   - *NaN poisoning*: a single skipped fold turned `meanFeatureWeight` into all-NaN in PLSDA and
     ENet, and ENet's `featureStability` counted NaN columns as "not selected".
   - *ENet used three different lambda paths* for the point estimate, the CI and the null; unified
     on `opts.lambdaGrid` via `enetLambdaGrid`.
   - Missing small-sample guards in the PLSDA quick-CV helpers; `quickCV_ENet` returning 0.5 rather
     than NaN on total failure; the PR permutation histogram plotting the ROC vector; a local
     variable named `diag` shadowing the builtin.
   - `opts.nrepeats` in both calling scripts was never read (the pipelines read `nRepeats`), so the
     option silently did nothing. `warnUnknownOptions` now catches this class of bug.

**Also fixed (separate commit): `group_tfce_from_subject_maps.m` permutation schemes.** Both designs
were producing invalid inference, and neither problem was confined to covariates.

- *One-sample with covariates could never reject.* The nuisance design carried a column of ones, but
  for a one-sample test the intercept **is** the effect of interest, so the group mean was swept into
  `Y_hat` and Freedman-Lane added it back into every permuted dataset. Covariates are now
  mean-centered and regressed out without an intercept column, via `residualizeFold`.
- *Two-sample had almost no power, with or without covariates.* `perm_idx` was applied to both the
  residuals and the group labels, which preserves their pairing and so permutes nothing: the permuted
  group sets are the observed ones, reordered. The t-statistic is then bit-identical to the observed
  (verified 200/200). End to end through pTFCE the maps differ only in the last bits, because
  floating-point summation is not associative — the null took **2 distinct values across 40
  permutations**, sd 0.011 against 0.608 for the fixed version, so p-values were decided by rounding
  noise. Labels are now held fixed while the residuals are permuted.
- `send(q,1)` ran unconditionally although `q` only existed when `'parallel'` was true, so that
  documented option errored instead of disabling parallelism. Now guarded and honoured via `parfor`'s
  worker cap.

Measured on synthetic data with a known answer (n=40, 200 datasets; nominal 0.05 has a 95% interval
of [0.020 0.080]):

| design | before → after, false positives | mean p | power |
|---|---|---|---|
| one-sample | 0.000 → **0.055** | 0.972 → 0.502 | 0.00 → **1.00** |
| two-sample | 0.060 → **0.040** | 0.475 → 0.532 | 0.10 → **1.00** |

**Any TFCE result from this function predates a working permutation test and should be re-run.**

Re-verified end to end against real CANlab / SPM12 / pTFCE once those dependencies were available
(`fmri_data_st` lives in the `canlab_single_trials` repo, which must be on the path). On data
containing a genuine effect, `p_global` pre-fix vs post-fix:

| design | pre-fix | post-fix |
|---|---|---|
| one-sample, no covariates | 0.0244 | 0.0244 |
| one-sample, **with covariates** | **1.0000** | **0.0244** |
| two-sample, no covariates | 0.4634 | **0.0244** |
| two-sample, with covariates | 0.4146 | **0.0244** |

The `send(q,1)` bug also reproduced live: the pre-fix code dies with `Unrecognized function or
variable 'q'` under `'parallel',false`.

*Constructing the null correctly is subtler than it looks, and cost two rounds of measurement here.*
The covariate-adjusted mean is only zero if the covariate effect is generated around the covariate's
own mean — `Y = b*(cv - mean(cv)) + e`. Validating against `Y = b*cv + e` makes every scheme look
broken, because the adjusted mean is then `b*mean(cv)`, nonzero in any finite sample.

**Results produced before this work are not numerically comparable to results produced after it**,
and prior ENet results in particular were selected by a flawed tuning rule.

### Audit of the remaining `functions/` files — findings, not yet fixed

The overhaul above covered the ML pipeline family and `group_tfce_from_subject_maps.m`. The
remaining files were read through afterwards. Nothing below is fixed yet.

**Silently wrong output — findings 1 and 2 are FIXED; finding 3 turned out to be something else.**

1. ~~**Atlas label gaps mis-map features to brain regions**~~ (all three plotters). **FIXED** via the
   shared `validateAtlasLabels`, which errors on labels missing from `1:p` and warns on labels beyond
   `p`. Verified end to end through the real plotter: a contiguous atlas runs clean, `[1 2 5 7 9]`
   with `p=5` is now blocked, and an atlas with extra regions runs with a warning. Original defect
   description follows. The painting loop is
   `for i = 1:p, mask = atlasData == i`, so the convention "atlas label `i` = column `i`" must hold
   exactly. The only guard is `if max(labels) < p, warning(...)`, which does **not** fire when the
   labels are non-contiguous. With labels `[1 2 5 7 9]` and `p = 5`: `max(labels)=9`, so no warning;
   features 3 and 4 paint onto nothing and are silently dropped, and feature 5 paints onto atlas
   label 5, which is the *third* ROI. Gaps arise whenever an ROI ends up with no voxels after
   masking or resampling. The check that catches it is `isequal(labels(:)', 1:p)`.
2. ~~**`thresholded_fmri_data_from_statistic_image` / `_from_pval_nii` threshold by value, not by
   significance.**~~ **FIXED** via the shared `maskToSignificant`, which parks non-significant voxels
   on a sentinel below the retained window so the range threshold selects exactly the significant
   set, and replaces the hardcoded `0.001` pad with a relative tolerance. Measured through CANlab's
   real `threshold` method on a case with 159 significant voxels: the old path kept **2544** voxels
   (2385 of them not significant), the new path keeps **exactly 159**. Original description follows. Both compute the significant voxels, then call
   `threshold([min(sig)-0.001 max(sig)+0.001], 'raw-between', ...)` — which keeps *every* voxel whose
   statistic falls in that range, significant or not. This is only equivalent to masking when the
   p-value is monotonic in the statistic, which it is not for voxelwise permutation p-values (each
   voxel has its own null). The `0.001` pad is also a hardcoded absolute tolerance on an arbitrary
   statistic scale.
3. ~~**`call_pTFCE` calls pTFCE with entirely the wrong arguments.**~~ **RESOLVED** by retiring
   pTFCE altogether (see "The TFCE stack" below). For the record, what was wrong: the installed
   signature is `pTFCE(imgZ, mask, Rd, V, Nh, Zest, C, verbose)` while `call_pTFCE` passed
   `(vol, voxsize, H, E, conn)`; pTFCE's `switch nargin` handles only cases 1-4, so the 6- and
   5-argument attempts always failed and every call fell through to the 4-argument form, where
   `mask` got a scalar voxel size (never used in pTFCE's body) and `Rd`/`V` got H = 2 and E = 0.5
   instead of a resel count and a voxel count. pTFCE has no `H` or `E` parameter at all, and `conn`
   never reached it. Measured impact: within a map the wrong parameterisation was exactly
   rank-preserving (Spearman 1.0000) so the spatial pattern survived, but values were inflated
   (max 7.00 vs 4.11) and across maps the rank correlation was 0.9654, so permutation p-values were
   approximately but not exactly right. **The same bug existed inline in
   `decoding_toolbox/LaBGAScore_decoding_SVM_between_subjects.m` and is fixed in the same change.**
   All TFCE numbers change; results predating this must be re-run.

**Hard errors — all now FIXED.**

0. `TopN` was not capped by the feature count in the PLSR and PLSDA plotters. Fixed; see the note
   above.
4. `tfce_one_fmri_dat`'s `switch sidedness` had no `otherwise`, leaving `tf` undefined. Resolved by
   the TFCE overhaul: that switch now lives in `tfce_volume`, which validates both `sidedness` and
   `tail` with identified errors.
5. `region_table.minP = []` errored when the table had no `minP` column, because assigning `[]`
   deletes a variable and deleting a non-existent one is an error rather than a no-op. Now guarded
   with `ismember('minP', ...)`.
6. Both `thresholded_fmri_data_from_*` preallocated with `height(region_table)` but looped over
   `size(region_obj,2)`. If those differed the array grew silently and the write-back failed with a
   table-height error. Both now use `numel(region_obj)`, and a genuine mismatch warns and returns
   the un-augmented table rather than killing a visualisation helper.
7. `thresholded_fmri_data_from_pval_nii` shadowed the builtin `path`. Fixed (renamed `out_dir`) as
   part of the significance-thresholding change.
8. No plotter checked `numel(roiNames) == p`. All three now error up front, before any NIfTI is
   written. Verified: 3 or 8 names against 5 features are both rejected.
9. `dice_statistic_image_by_roi` checked nothing. It now verifies the input types, that the two
   statistic images share a voxel space, that `roi_img` is in that same space, that supplied
   `binary_mask` entries match the voxel count, and that p-values exist when no mask is given; and
   its option `switch` gained the `otherwise` error its sibling always had, so a mistyped option no
   longer falls through to the default. Its accumulators are also preallocated rather than grown.

**Dead code — resolved.** `pvals_from_ranks` has been deleted: it had no callers in any repo on
this machine, it assigned arbitrary distinct p-values to tied statistics (a real problem for TFCE
maps, which are full of ties at zero), and despite its comment it computed a *spatial rank across
voxels* rather than a p-value against a null.

The two `dice_*` functions are kept and corrected. They had drifted apart: `dice_statistic_image`
thresholded with `p <= threshold` while `dice_statistic_image_by_roi` used `p < threshold`, so the
pair disagreed on voxels sitting exactly at the threshold. Both now use strict `<`, matching
`thresholded_fmri_data_from_*`. `dice_statistic_image` also lost a commented-out "repair geometry"
block and gained the binary-mask size check its sibling has; `dice_statistic_image_by_roi` gained
guards against being handed a continuous map instead of a label image. Verified that the whole-image
and per-ROI counts agree over the labelled voxels.

Note on earlier claims in this file: "no callers anywhere in the repo" was twice scoped to
LaBGAScore only. `tfce_transform_3d` is called by
`CanlabCore/Statistics_tools/searchlight_disti_Lukas.m`, and
`thresholded_fmri_data_from_statistic_image` by
`CANlab_help_examples/.../c2_SVM_contrasts_masked.m`. Check every repo before calling something dead.


## Three layers, different contracts

- `scripts/` — **templates meant to be copied into a study's `code` subdataset and edited**, not
  library code. They are not standalone: they `load` `.mat` files produced by an earlier pipeline
  step and assume workspace variables set by other scripts (see "Cross-repo data flow"). Editing a
  script here changes the template for future studies, not any running analysis.
- `functions/` — generic, study-agnostic library code, called from `scripts/` here and from adapted
  copies in study repos. Safe to fix in place; behavior changes propagate to every caller.
- `classes/ProgressTracker.m` — the repo's only class. A `handle` progress/ETA printer that uses
  `drawnow limitrate` so it prints from inside `parfor`. Used by
  `functions/group_tfce_from_subject_maps.m` and by `decoding_toolbox/` one level up.

## The ML pipeline family

Four pipelines and three plotters form the largest subsystem here:

| Pipeline | Outcome | Estimator | Plotter |
|---|---|---|---|
| `PLSR_neuroimaging_pipeline` | continuous | `plsregress` | `plot_PLSR_diagnostics_neuroimaging` |
| `PLSDA_neuroimaging_pipeline` | binary | `plsregress` | `plot_PLSDA_diagnostics_neuroimaging` |
| `PLSDA_paired_neuroimaging_pipeline` | binary, within-subject | `plsregress` | (reuses the PLSDA plotter) |
| `ENet_neuroimaging_pipeline` | binary | `lassoglm` | `plot_ENet_diagnostics_neuroimaging` |

All take a `[n x p]` subjects × features matrix `X` (PET ROI binding, fMRI ROI betas, morphometry,
connectivity edges, graph metrics) plus an `opts` struct, and return one `results` struct.

**Shared skeleton, now genuinely shared.** Each pipeline body runs the same numbered stages —
`0. Defaults` → `1. Outcome preparation` → `2. Repeated nested CV` → global-signal baseline → final
model → stability metrics → `Permutation testing` → `Bootstrap CI (out-of-bag)` → `Learning curve`.
The per-fold helpers used to be local subfunctions copied into each file (`applyScaling` existed in
six copies, `capLV` in three). They are now **standalone files in `functions/`**, so a fix reaches
every caller:

| helper | role |
|---|---|
| `foldPreprocess` | the one entry point every fold uses: residualize, then scale |
| `residualizeFold` / `residualizeY` | train-only nuisance regression (truncated SVD, rank-safe) |
| `applyScaling` | train-only z-score / center / none |
| `capLV` | latent-variable cap; call on the **post**-preprocessing matrix |
| `validateCovariates` / `warnUnknownOptions` | up-front input and option checks |
| `quickCV_PLSR` / `quickCV_PLSDA{,_PR}` / `quickCV_PLSDA_core` / `quickCV_ENet{,_PR}` | fast tuned CV for permutation, learning curve, matched observed |
| `bootstrapOOB_{PLSR,PLSDA,ENet}` | one out-of-bag bootstrap replicate |
| `selectENetHyperparams` / `enetLambdaGrid` | fold-averaged ENet tuning on one shared lambda path |
| `globalBaselineCV` | cross-validated companion to the in-sample global baseline |
| `setParforStream` | per-iteration substreams so `parfor` results are reproducible |
| `makeGroupedFolds` / `swapWithinSubjectLabels` / `quickGroupedCV` | paired-design subject blocking |

Note the naming: `quickCV` could not simply be lifted out, because PLSR's returned Q² and PLSDA's
returned AUC under the same name. They are split by estimator, matching the pre-existing
`quickCV_ENet` convention.

**Invariants the whole family relies on:**

- *Leakage-free preprocessing* — every fold calls `foldPreprocess`, which fits both the nuisance
  coefficients and the scaling constants on the training fold alone. Any refactor that centers,
  z-scores, or residualizes `X` before the CV loop silently invalidates every performance estimate
  downstream. This applies at **five** sites per pipeline: outer fold, inner tuning fold,
  `quickCV_*`, `bootstrapOOB_*`, and the interpretation-only final model.
- *Covariates travel as an explicit argument, never inside `opts`* — helpers take `C` alongside `X`.
  Bootstrap resampling and learning-curve subsampling reindex subjects, and a full-length covariate
  matrix reaching a subsetted `X` would misalign silently instead of erroring. `foldPreprocess`
  asserts matching row counts to turn that class of bug into a hard failure.
- *Order is residualize-then-scale* — decided once, in `foldPreprocess`, and nowhere else.
  Residualization lowers `rank(X)`, so `capLV` must see the post-residualization matrix.
- *Positive class = `max(Y)`* — binary pipelines binarize as `double(Y == max(Y))`.
- *Feature order is the atlas index* — atlas integer labels `1..p` must correspond to the column
  order of `X`. This convention is load-bearing across pipeline → plotter → atlas file; the plotters
  paint `results` values into a NIfTI by treating column `i` as atlas label `i`.
- *`results` is the interface* — the plotters are decoupled from the pipelines and consume only
  documented `results` fields (`VIP`, `stabilityZ`, `finalXLoadings`, selection frequencies, …),
  degrading gracefully when optional ones are absent. Adding a field is safe; renaming one breaks
  the plotter silently. The paired PLS-DA pipeline exists to feed the PLSDA plotter the same field
  names from a subject-blocked CV scheme (each subject appears exactly twice, once per class — it
  errors otherwise).

## The TFCE stack

Layered, with an external dependency at the bottom:

```
group_tfce_from_subject_maps   % group inference: sign-flip / label-exchange permutation,
                               % optional Freedman-Lane covariate control, parfor over permutations
  └─ tfce_one_fmri_dat         % parfor-safe single-map transform; rebuilds a 3D volume from the
                               % masked vector and hands it to tfce_volume
       └─ tfce_volume          % sidedness/tail handling and the shared integration step
            └─ tfce_transform_3d   % classic TFCE (Smith & Nichols 2009), via bwconncomp
```

`decoding_toolbox/LaBGAScore_decoding_SVM_between_subjects.m` enters the same stack at
`tfce_volume`, so there is one TFCE implementation in the repo rather than two.
`pvals_from_ranks.m` turns permutation ranks into uncorrected voxelwise p-values.

**One algorithm, deliberately.** The repo previously mixed classic TFCE with pTFCE, and the pTFCE
side was calling its dependency with mismatched arguments. Standardising on classic TFCE was chosen
because: `CanlabCore/Statistics_tools/searchlight_disti_Lukas.m` already uses `tfce_transform_3d`
correctly for the searchlight path; the options the whole ecosystem exposes (`H`, `E`, `conn`) are
classic-TFCE parameters that pTFCE does not have; pTFCE exists to avoid permutation, so using it
inside a permutation loop discards its only advantage while paying its costs (Z-map, resel count,
GRF assumptions); and it is 10-13x faster in that loop (measured 0.50 s vs 6.36 s per map on a
64x76x64 volume, i.e. 8 min vs 106 min at nPerm = 1000). `call_pTFCE.m` has been deleted and no
`pTFCE` call remains in the repo.

**The integration step is shared.** `tfce_transform_3d` defaults to `dh = max(stat)/100`, a grid
derived from each map's own maximum, which inside a permutation loop gives every permuted map a
different grid. The observed map now fixes `dh` and all permutations reuse it (recorded as
`info.tfce_dh`). Measured effect of not doing this: ~0.6% difference in the TFCE maximum and a rank
correlation of 0.995, small but injected straight into the null.

**Turning p-maps into something viewable:** `thresholded_fmri_data_from_statistic_image.m` (from an
in-memory `statistic_image`) and `thresholded_fmri_data_from_pval_nii.m` (from a p-value `.nii` on
disk) both take the p-map plus the *original statistic vector*, and return a thresholded
`fmri_data_st`, a `region` object, and a labeled table. They exist because CANlab's
`statistic_image.threshold()` does not handle non-standard statistics (AUC, TFCE); both do the
masking manually. Both return `[]` for all outputs when nothing survives threshold — callers must
check. `dice_statistic_image.m` / `dice_statistic_image_by_roi.m` compare two thresholded maps.

## Cross-repo data flow

`scripts/` here sit **downstream of the `CANlab_help_examples` (LaBGAS fork) second-level templates**
and read their outputs:

- `LaBGAScore_secondlevel_roi_run_plot_PLS_ENet_pipeline.m` loads `image_names_and_setup.mat`,
  `roi_stats_<mygroupnamefield>_<scaling_string>_<results_suffix>.mat` (from
  `prep_3a_run_second_level_regression_and_save.m`) and the combined-ROI `.mat` from
  `LaBGAScore_atlas_rois_from_atlas.m`, then calls the pipelines above.
- `LaBGAScore_secondlevel_extractparcels_sessions.m` likewise loads `prep_3a` results.

Consequence: option variables such as `mygroupnamefield`, `results_suffix`, `myscaling_glm`,
`atlasname_glm`, `roi_modelname`, and `roi_set_name` **must match the values used in the upstream
`prep_3a` / `a2_set_default_options` run**, or the `load` finds no file or the wrong one. Directories
come from `prep/LaBGAScore_prep_s0_define_directories.m`. Two helpers used by these scripts live
outside this folder: `clean/LaBGAScore_smart_parallel_pool_setup.m` and
`figures/save_all_open_figures_smart.m`.

The remaining scripts stand apart: `LaBGAScore_secondlevel_MS_mat_pipeline.m` builds an SPM batch for
the MACS toolbox model space; `LaBGAScore_secondlevel_mvpa_beta_maps_conn.m` runs MVPA regression on
CONN-toolbox connectivity maps (explicitly *not* yet adapted to the standard LaBGAS dataset layout);
`LaBGAScore_secondlevel_ooFmriDataObjML_example.m` is illustrative only and assumes `fmri_dat` and
several options already exist in the workspace.

## Dependencies used in this folder

Beyond CANlab (`fmri_data_st`, `statistic_image`, `atlas`, `region`) and SPM12: **Statistics and
Machine Learning Toolbox** (`plsregress`, `lassoglm`, `perfcurve`, `cvpartition`, `fitglm`,
`robustfit`) is required by every ML pipeline; **Image Processing Toolbox** (`bwconncomp`) by
`tfce_transform_3d`; **Parallel Computing Toolbox** by the `parfor` loops in the three main pipelines
and the TFCE group function. The plotters call SPM directly (`spm_vol`, `spm_read_vols`,
`spm_write_vol`) to write NIfTI maps. External, not vendored: **pTFCE**, **MACS toolbox**,
**ooFmriDataObjML**, **CONN**.

## Checks

No build system, test suite, or CI. Static analysis over this subtree only:

```bash
matlab -batch "addpath(genpath('..')); LaBGAScore_check_all_scripts(pwd)"
```

`checkcode` catches parse errors but not undefined variables, missing functions, or wrong indexing —
the defect classes that actually occur here. Read the code path; don't rely on it.

## Conventions

Files in `functions/` use MATLAB's standard function help-text convention (H1 line, `USAGE`,
`INPUTS`, `OUTPUTS`, `See also`). Files in `scripts/` use the repo's script-header convention
(`%% scriptname.m` title, asterisked `*USAGE*` / `*OPTIONS*` / `*DEPENDENCIES*` / `*NOTES*` sections,
dashed separators, then an author/date/version block ending in a `last modified: YYYY/MM/DD` line).
Script `*OPTIONS*` blocks document the user-editable variables in the `USER SETTINGS — EDIT THESE`
section at the top of the body; keep the two in sync when adding an option.
