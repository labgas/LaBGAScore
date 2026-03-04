# plot_PLSDA_diagnostics_neuroimaging — User Guide

`plot_PLSDA_diagnostics_neuroimaging` is a generic visualization/diagnostics helper for **PLS-DA** results produced by the companion pipeline (`PLSDA_neuroimaging_pipeline` / `TSPO_PLSDA_pipeline`).

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
    'LV', 2, 'OutPrefix', 'PLSDA_run1');

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
- `results.selectionFrequency` **[p×1]** — top-K selection frequency (for top-K tables)
- `results.featureStability` **[p×1]** — proportion non-zero across runs (stem plot)

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
| `TopN` | 20 | Number of rows to display at top of ROI table (table file exports all rows). |
| `TopK` | 20 | Number of top features (by absolute meanFeatureWeight) for bar/table export. |
| `LV` | 2 | LV index to visualize as a brain map. |
| `OutPrefix` | `'PLSDA'` | Prefix for exported files (CSV + NIfTI). |
| `VIP_thresh` | 1 | VIP threshold for “robust contributor” marking/labeling. |
| `stab_thresh` | 2 | stabilityZ threshold for “robust contributor” marking/labeling. |
| `MapPrctile` | 70 | Threshold percentile for LV map (keeps top `100-MapPrctile`% by |loading|). |
| `RelaxIfEmpty` | true | If no ROI passes robust thresholds, relax thresholds (75th percentile) for visualization only. |

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

4. **Feature stability stem plot** (optional)
Requires `results.featureStability`.

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
- When p is large and n small, stabilityZ may be noisy; prefer robust thresholds plus selectionFrequency/signStability outputs from the pipeline.
- When p is moderate and n larger, stabilityZ becomes more interpretable; you may tighten thresholds.

---

## Minimal reproducible call

```matlab
ROI_table = plot_PLSDA_diagnostics_neuroimaging(results, results.finalXLoadings, roiNames, atlasFile, ...
    'LV', 2, 'OutPrefix', 'Study1_PLSDA', 'TopN', 20, 'TopK', 20);
```

