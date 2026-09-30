# secondlevel — group-level statistics, MVPA/ML pipelines, TFCE, and reporting

The library the group-level analyses are built from: 48 functions, 7 template
scripts, and the repo's only class. Most of it is called from the second-level
templates in the [`CANlab_help_examples`](https://github.com/labgas/CANlab_help_examples)
LaBGAS fork rather than run directly from here.

This is an **index**. Seven standalone guides already cover the ML pipelines in
depth, and `CLAUDE.md` in this folder covers the design and audit history — both
are authoritative for their subjects and are not restated here.

---

## Where to look for what

| you want | read |
|---|---|
| to run PLSR, PLS-DA or Elastic Net on neuroimaging data | the pipeline guide for that method (table below) |
| to understand or adapt the diagnostic plots | the plotting guide for that method |
| the design of the ML family, the layer contracts, the audit history | [`CLAUDE.md`](CLAUDE.md) in this folder |
| how a study's second level is actually driven | that repo's `Second_level_analysis_template_scripts/README.md` |
| everything else in here | this file |

**The seven ML guides:**

| method | pipeline | plotting |
|---|---|---|
| Elastic Net | [`README_ENet_neuroimaging_pipeline.md`](README_ENet_neuroimaging_pipeline.md) | [`README_ENet_plotting.md`](README_ENet_plotting.md) |
| PLS-DA | [`README_PLSDA_neuroimaging_pipeline.md`](README_PLSDA_neuroimaging_pipeline.md) | [`README_PLSDA_plotting.md`](README_PLSDA_plotting.md) |
| PLS-DA, paired | [`README_PLSDA_paired_neuroimaging_pipeline.md`](README_PLSDA_paired_neuroimaging_pipeline.md) | — |
| PLSR | [`README_PLSR_neuroimaging_pipeline.md`](README_PLSR_neuroimaging_pipeline.md) | [`README_PLSR_plotting.md`](README_PLSR_plotting.md) |

## Three layers, three different contracts

Spelled out in full in [`CLAUDE.md`](CLAUDE.md); the short version, because
mixing them up is the easiest mistake here:

- **`scripts/`** (7) — templates to be **copied into a study's `code` subdataset
  and edited**, not library code. They are not standalone: they `load` `.mat`
  files from an earlier pipeline step and assume workspace variables other
  scripts set. Editing one changes the template for future studies, not any
  running analysis.
- **`functions/`** (48) — generic, study-agnostic library code, called from
  `scripts/` here and from adapted copies in study repos. Safe to fix in place —
  and a behaviour change propagates to every caller, in every repo.
- **`classes/ProgressTracker.m`** (1) — the repo's only class. A `handle`
  progress/ETA printer that uses `drawnow limitrate` so it prints from inside
  `parfor`. Used by `group_tfce_from_subject_maps` and by `decoding_toolbox/`.

## What is in `functions/`

### The ML pipeline family (~30 files)

Four pipelines — `PLSR_neuroimaging_pipeline`, `PLSDA_neuroimaging_pipeline`,
`PLSDA_paired_neuroimaging_pipeline`, `ENet_neuroimaging_pipeline` — three
matching `plot_*_diagnostics_neuroimaging` plotters, and the shared machinery
they sit on: cross-validation (`quickCV_*`, `quickGroupedCV`, `makeGroupedFolds`,
`globalBaselineCV`), out-of-bag bootstrap (`bootstrapOOB_*`), in-fold
preprocessing (`foldPreprocess`, `applyScaling`, `residualizeFold`,
`residualizeY`), hyperparameter selection (`selectENetHyperparams`,
`enetLambdaGrid`, `capLV`), and input validation (`validateCovariates`,
`validateAtlasLabels`, `warnUnknownOptions`).

Read the guides above before touching any of these. The one invariant worth
stating here: **preprocessing that could see the outcome happens inside the
fold**, which is what `foldPreprocess` and `residualizeFold` exist to enforce.

### The TFCE stack (5 files)

`tfce_transform_3d` → `tfce_volume` → `tfce_one_fmri_dat`, with
`group_tfce_from_subject_maps` for group TFCE from subject-level maps and
`tfce_fwe_from_null` for max-statistic FWE p-values out of a saved permutation
null. `tfce_volume` is **shared with `decoding_toolbox/`**, so its inference and
the SVM pipeline's are the same code.

TFCE already integrates cluster extent into the statistic, so **do not add a
cluster-extent threshold on top of it** — that removes signal the correction has
already accounted for.

### Thresholding, tables and display (7 files)

| function | |
|---|---|
| `thresholded_fmri_data_from_statistic_image.m` | thresholded `fmri_data` from a `statistic_image` |
| `thresholded_fmri_data_from_pval_nii.m` | the same from a p-value NIfTI |
| `LaBGAScore_region_table.m` | `@region/table` vendored with two upstream defects removed — see below |
| `LaBGAScore_region_table_safe.m` | three-output wrapper, so one undisplayable contrast does not abort a whole report |
| `LaBGAScore_blob_montage.m` | overview + regioncenters montages. Display only: it draws, it does not threshold, label or tabulate. Skips regioncenters above 21 regions, where one titled panel per region stops being readable |
| `dice_statistic_image.m`, `dice_statistic_image_by_roi.m` | Dice overlap between two thresholded maps, whole-image and per-ROI. Both threshold with strict `<` and return **`NaN`, not 0**, where neither map has a suprathreshold voxel |

**Two things about `LaBGAScore_region_table` that change what a number means:**

1. Upstream `@region/table` ends `get_signed_max` with `maxZ = norminv(1 - 1E-12)`
   — **7.0345** — and clips every value above it, not only the infinities its own
   comment describes. This copy does not.
2. It relabels the Bayes column. `@statistic_image/estimateBayesFactor` returns
   **2·ln(BF10)**, not BF, so a column headed `maxBF` showing 15.05 was a Bayes
   factor of 1854. The column is now `max_2lnBF`, with `maxBF10 = exp(max_2lnBF/2)`
   beside it.

Those two compounded: clipping at 7.0345 on a 2·ln scale caps reported evidence
at BF10 ≈ 34, so an unclipped peak of 39591 printed as 7.03. **Converting a
stored value: `BF10 = exp(v/2)`, never `exp(v)`.**

### PDM (3 files)

`LaBGAScore_pdm_report` and `LaBGAScore_pdm_report_image` in `functions/`, plus
`LaBGAScore_pdm_regenerate_reports.m` at the top level for rebuilding reports
from existing results.

## `scripts/`

| script | |
|---|---|
| `LaBGAScore_secondlevel_extractparcels_sessions.m` | extracts parcel-wise values from a parcellation for mixed-model analysis |
| `LaBGAScore_secondlevel_extractclusters_sessions.m` | the voxel-based counterpart, working on `c2a` output (FDR, TFCE FDR or TFCE FWE) |
| `LaBGAScore_secondlevel_roi_run_plot_PLS_ENet_pipeline.m`, `..._PLSR_pipeline.m` | ROI-level runners for the ML pipelines |
| `LaBGAScore_secondlevel_MS_mat_pipeline.m` | multi-session `.mat` pipeline |
| `LaBGAScore_secondlevel_mvpa_beta_maps_conn.m` | MVPA beta maps from CONN output |
| `LaBGAScore_secondlevel_ooFmriDataObjML_example.m` | worked example of CANlab's object-oriented ML interface |

## Cross-repo data flow

Five functions here are called from the `CANlab_help_examples` LaBGAS fork —
`c2a_second_level_regression`, `prep_3a_run_second_level_regression_and_save`
and `c2_SVM_contrasts_masked` between them use `group_tfce_from_subject_maps`,
`thresholded_fmri_data_from_statistic_image`, `tfce_fwe_from_null`,
`LaBGAScore_region_table` and `LaBGAScore_region_table_safe`. `tfce_volume` is
additionally called from `decoding_toolbox/` in this repo. The root
[`DEPENDENCIES.md`](../DEPENDENCIES.md) is the generated, authoritative record,
and the cross-repo table in the root [`README.md`](../README.md) is the readable
version.

Two consequences: **a behaviour change in `functions/` reaches study analyses
without anything here being re-run**, and "no callers in this repo" is not
"dead" — check every clone before deleting anything (`tfce_transform_3d` is
called from CanlabCore, and `thresholded_fmri_data_from_statistic_image` from
the fork).

## Dependencies

**CanlabCore** throughout (`fmri_data`, `statistic_image`, `region`, `atlas`,
`canlab_results_fmridisplay`), **Neuroimaging_Pattern_Masks** for atlases,
**SPM12** for image I/O, and MATLAB's **Statistics and Machine Learning
Toolbox** for the ML family. Neither CANlab repo is vendored.

## Status

The ML pipeline family and `group_tfce_from_subject_maps` were overhauled, and
the remaining `functions/` files were then audited — nine hard errors and three
silently-wrong-output findings, **all now fixed or resolved**, with
`secondlevel/functions/` reporting no parse errors and no unsuppressed Code
Analyzer messages. [`CLAUDE.md`](CLAUDE.md) has the itemised history.
