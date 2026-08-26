# plot_PLSDA_diagnostics_neuroimaging — User Guide

`plot_PLSDA_diagnostics_neuroimaging` is a generic visualization/diagnostics helper for **PLS-DA** results produced by the companion pipeline `PLSDA_neuroimaging_pipeline` (or a study-specific adaptation of it).

It creates:
- **tables** summarizing ROI importance and robustness
- **scatter plots** (VIP vs stabilityZ)
- **NIfTI maps** for VIP, stabilityZ, and LV loadings
- **multi-slice figures** for quick neuroanatomical inspection
- optional **top-K weight tables/plots** and **stability stem plots** (if present in `results`)

This function is designed to be copy-paste friendly for neuroimaging ROI/parcel data where an atlas uses integer labels `1..p` matching the feature order in `X`.

---

## Quick start

```matlab
% After running the pipeline:
% results = PLSDA_neuroimaging_pipeline(X, Y, opts);

% Use final model loadings if you have them:
ROI_table = plot_PLSDA_diagnostics_neuroimaging(results, results.finalXLoadings, roiNames, atlasFile, ...
    'LV', 1, 'OutPrefix', 'PLSDA_run1');

% Or provide XL explicitly (e.g., from your own plsregress call):
ROI_table = plot_PLSDA_diagnostics_neuroimaging(results, XL, roiNames, atlasFile);
```

---

## Inputs

### `results` (struct)
PLS-DA results struct, expected fields:
- `results.VIP` **[p×1]** — variable-importance-in-projection scores from the final model
- `results.stabilityZ` **[p×1]** — stability statistic from CV weights (meanBeta / sdBeta)

Optional (used if present):
- `results.finalXLoadings` **[p×LV]** — final model X-loadings (`XL`)
- `results.finalLV` — final number of LVs (used to clamp requested `LV`)
- `results.meanFeatureWeight` **[p×1]** — mean CV beta weights (for top-K plots/tables)
- `results.signStability` **[p×1]** — proportion of runs matching mean sign (top-K table column, stem plot)

### `XL` (numeric, [p×LV])
X-loadings matrix from `plsregress`, used to map LV loadings to ROIs.
- If `XL` is empty, the function tries `results.finalXLoadings`.

**Important:** `XL` must correspond to the same feature ordering as your columns in `X`.

### `roiNames` (cellstr/string, length p)
ROI/feature labels. If empty, defaults to `ROI_001 ... ROI_p`.

### `atlasFile` (char/string)
Path to a labeled atlas NIfTI used for ROI mapping.
- Atlas should contain integer labels `1..p` (0 as background).
- The function maps `feature i` → all voxels where `atlasData == i`.

---

## Name–Value Options

| Option | Default | Meaning |
|---|---:|---|
| `TopN` | 20 | Number of rows to display at top of ROI table, and number of top features (by absolute meanFeatureWeight) shown in the bar/table export — this single option controls both (table file exports all ROI rows regardless). |
| `LV` | 1 | LV index to visualize as a brain map. |
| `OutPrefix` | `'PLSDA'` | Prefix for exported files (CSV + NIfTI). |
| `VIP_thresh` | 1 | VIP threshold for “robust contributor” marking/labeling. |
| `stab_thresh` | 2 | stabilityZ threshold for “robust contributor” marking/labeling. |
| `MapPrctile` | 70 | Threshold percentile for LV map (keeps top `100-MapPrctile`% by |loading|). |
| `RelaxIfEmpty` | true | If no ROI passes robust thresholds, relax thresholds (75th percentile) for visualization only. |
| `UnderlayFile` | `''` | Optional structural underlay NIfTI for the multi-slice figure. Must match `atlasFile`'s voxel space/dimensions. |

---

## Outputs

### `ROI_table` (MATLAB table)
Returned table with columns:
- `ROI` — ROI/feature label
- `VIP` — VIP score
- `stabilityZ` — stabilityZ statistic
- `RobustContributor` — boolean, based on thresholds

It is also saved to disk as:
- `'<OutPrefix>_ROI_VIP_stability.csv'`

---

## Files created (exports)

1. **CSV**
- `'<OutPrefix>_ROI_VIP_stability.csv'`  
  Full ROI table (VIP, stabilityZ, robust flag).
- `'<OutPrefix>_top<k>_weights.csv'` (if `results.meanFeatureWeight` exists)  
  Top-K weights table.

2. **NIfTI maps** (requires SPM)
- `'<OutPrefix>_VIP_map.nii'` — VIP values mapped to atlas
- `'<OutPrefix>_stabilityZ_map.nii'` — stabilityZ mapped to atlas
- `'<OutPrefix>_LV<LV>_map.nii'` — LV loadings mapped to atlas (normalized to max |loading|)
- `'<OutPrefix>_LV<LV>_map_thresh.nii'` — thresholded LV map based on `MapPrctile`

---

## Figures created (interactive MATLAB figures)

1. **Scatter plot: VIP vs stabilityZ**
- X-axis: VIP
- Y-axis: stabilityZ
- Horizontal/vertical lines at `VIP_thresh` and `stab_thresh`
- Robust ROIs labeled with text.

2. **Multi-slice 3×N panel**
Rows:
- LV map (selected LV)
- VIP map (thresholded at `VIP_thresh`)
- stabilityZ map (thresholded at `stab_thresh`)
Robust ROI labels are overlaid on each slice.

3. **Top-K bar plot** (optional)
Requires `results.meanFeatureWeight`:
- Bars show mean weight for the top-K absolute weights.

4. **Sign stability stem plot** (optional)
Requires `results.signStability`.

---

## Robust contributor definition

By default:
- `RobustContributor = (VIP > VIP_thresh) & (abs(stabilityZ) > stab_thresh)`

This mirrors a common PLS-DA reporting convention:
- VIP captures *importance in the final fitted model*
- stabilityZ captures *robustness across CV re-fits*

### When `RelaxIfEmpty = true`
If **no** ROI passes the thresholds, the function relaxes them **for visualization only**:
- `VIP_thresh_vis = prctile(VIP, 75)`
- `stab_thresh_vis = prctile(abs(stabilityZ), 75)`

The exported table still contains the robust flag computed with the *visualization thresholds* in this scenario (because `isRobust` is updated).
If you want strict thresholds no matter what, set `RelaxIfEmpty = false`.

---

## Atlas/feature alignment checklist

This function assumes:
1. The **column order** of your feature matrix `X` matches:
2. The **row order** of `results.VIP`, `results.stabilityZ`, and `XL`, and:
3. The atlas labels `1..p` correspond to those same features.

Common failure modes:
- Atlas uses non-contiguous labels (e.g., 1001..2000) → mapping will be wrong.
- Atlas has fewer labels than `p` → function warns if `max(atlasLabel) < p`.
- ROI list includes background/0 index → should be removed upstream.

If your atlas labels are not 1..p, remap them to 1..p before using this function.

---

## Interpreting the LV map

- LV loadings represent a **multivariate pattern** (a component) rather than independent univariate effects.
- Positive and negative loadings indicate opposite sides of the LV contrast.
- If your positive class is coded as `Y == max(Y)`, the **direction** of the component relative to classes depends on the model/sign conventions; treat LV maps as *pattern descriptors*, not direct “increase/decrease” claims without verifying sign using class means on scores.

---

## Dependencies

Requires:
- **SPM** (for atlas I/O): `spm_vol`, `spm_read_vols`, `spm_write_vol`
- MATLAB Statistics and Machine Learning Toolbox (already needed by pipeline): tables/plotting basics

---

## Recommended use by study type

### ROI-level features (PET binding, fMRI betas, morphometry)
- Atlas mapping is usually straightforward.
- `VIP_thresh=1` and `stab_thresh=2` are reasonable starting points.

### Connectivity edges / graph metrics
- If your “features” are edges, you **cannot** map directly to an atlas unless you have an edge-to-voxel strategy.
- You can still use:
  - ROI table
  - scatter plot
  - top-K plots
  but skip atlas/NIfTI steps by not calling this function or by adapting it to your visualization needs.

### p/n ratio considerations
- When p is large and n small, stabilityZ may be noisy; prefer robust thresholds plus the pipeline's signStability output.
- When p is moderate and n larger, stabilityZ becomes more interpretable; you may tighten thresholds.

---

## Minimal reproducible call

```matlab
ROI_table = plot_PLSDA_diagnostics_neuroimaging(results, results.finalXLoadings, roiNames, atlasFile, ...
    'LV', 1, 'OutPrefix', 'Study1_PLSDA', 'TopN', 20);
```

---

## New `results` fields (pipeline overhaul)

The pipelines now add the fields below. This plotter does not require any of
them and is unaffected if they are absent, but they are worth knowing when
reading a saved `results` struct:

| field | meaning |
|---|---|
| `covariateInfo` | struct: `.used`, `.names`, `.nCov`, `.rank`, `.residualizeY`, `.order`, `.permScheme`. Records whether fold-wise covariate regression was applied and how. |
| `quickCV_observed` | observed statistic computed with the same quick-CV estimator that generates the permutation null. `permutation_p` is now tested against this, not against the headline metric. |
| `seed` | the RNG seed actually used; results are reproducible from it and independent of parallel pool size. |
| `*_global_cv` | cross-validated companion to the in-sample global-signal baseline. Compare the model against the `_cv` value, not the in-sample one. |

No existing field was renamed or changed in meaning, so previously saved
`results` structs still plot correctly.

**Two caveats when comparing old and new outputs.** Permutation p-values from
before this change are inflated (the null used a different estimator than the
observed statistic, giving a measured false-positive rate near 40 % at
alpha = 0.05). And Elastic Net hyperparameters were previously selected by the
single best-performing inner fold rather than a fold-averaged score, which
produced markedly under-regularized, unstable weight maps.

---

## Input contracts now enforced

Three requirements that used to fail silently, or fail late, are checked up
front. All three abort **before** any NIfTI or figure is written.

### Atlas labels must be exactly `1..p`

The painting loop is `for i = 1:p, mask = atlasData == i`, so atlas integer
label `i` must BE feature `i`. `validateAtlasLabels` now errors when a label in
`1:p` is missing from the atlas volume, and warns when the atlas carries labels
beyond `p` (those regions simply stay empty in the maps).

The old guard only tested `max(labels) < p`, which misses the case that
actually bites: **gaps**. With labels `[1 2 5 7 9]` and `p = 5`, `max(labels)`
is 9 so nothing fired, yet features 3 and 4 matched no voxels and were dropped,
and feature 5 was painted onto atlas label 5 — the *third* ROI in the set. The
result looked like a perfectly normal brain map. Gaps arise whenever an ROI
loses all of its voxels during masking or resampling.

If you hit this error, rebuild the ROI atlas so its labels run `1..p` without
gaps, or subset `X` to match.

### `roiNames` must have exactly `p` entries

Too few names used to throw part way through, after the NIfTIs had been
written; too many, or names from a different ROI set, silently mislabelled
every row of `ROI_table`.

### `TopN` is capped at `p`

`TopN` is a display and export count. With the default `TopN = 20` and fewer
than 20 features the top-weights block used to throw "Index exceeds the number
of array elements", again after the NIfTIs were written. The exported CSV is
named for the true row count, so a 5-feature analysis writes
`..._top5_weights.csv`.
