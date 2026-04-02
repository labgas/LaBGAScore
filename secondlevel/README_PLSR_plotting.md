# plot_PLSR_diagnostics_neuroimaging — User Guide

`plot_PLSR_diagnostics_neuroimaging` is a generic visualization/diagnostics helper for **PLSR** results produced by the companion pipeline (`PLSR_neuroimaging_pipeline`).

It creates:

- **tables** summarizing ROI importance and robustness  
- **scatter plots** (VIP vs stabilityZ)  
- **NIfTI maps** for VIP, stabilityZ, and LV loadings  
- **multi-slice figures** for quick neuroanatomical inspection  
- optional **top-K weight tables/plots** and **stability stem plots** (if present in `results`)  
- optional **predicted vs observed plots** based on held-out CV predictions  
- optional **permutation**, **bootstrap**, and **learning curve** diagnostics (if present in `results`)

This function is designed to be copy-paste friendly for neuroimaging ROI/parcel data where an atlas uses integer labels `1..p` matching the feature order in `X`.

---

# Quick start

```matlab
% After running the pipeline:
% results = PLSR_neuroimaging_pipeline(X, Y, opts);

% Use final model loadings if you have them:
ROI_table = plot_PLSR_diagnostics_neuroimaging(results, results.finalXLoadings, roiNames, atlasFile, ...
    'LV', 1, 'OutPrefix', 'PLSR_run1');

% Or provide XL explicitly (e.g., from your own plsregress call):
ROI_table = plot_PLSR_diagnostics_neuroimaging(results, XL, roiNames, atlasFile);
```

---

# Inputs

## `results` (struct)

PLSR results struct, expected fields:

- `results.VIP` **[p×1]** — variable-importance-in-projection scores from the final model  
- `results.stabilityZ` **[p×1]** — stability statistic from CV weights (meanBeta / sdBeta)

Optional (used if present):

- `results.finalXLoadings` **[p×LV]** — final model X-loadings (`XL`)  
- `results.finalLV` — final number of LVs (used to clamp requested `LV`)  
- `results.meanFeatureWeight` **[p×1]** — mean CV beta weights (for top-K plots/tables)  
- `results.meanBeta` **[p×1]** — mean CV beta weights (preferred if present)  
- `results.selectionFrequency` **[p×1]** — top-K selection frequency (for top-K tables)  
- `results.signStability` **[p×1]** — sign consistency across CV runs  
- `results.featureStability` **[p×1]** — proportion non-zero across runs (stem plot)

Optional predictive diagnostics (from the pipeline):

- `results.cvObserved` **[N×1]** — observed Y values from outer CV  
- `results.cvPredicted` **[N×1]** — predicted Y values from outer CV  
- `results.cvRepeatID` **[N×1]** — repeat index for each held-out prediction  
- `results.cvSubjectID` **[N×1]** — subject index for each held-out prediction  

Optional resampling diagnostics:

- `results.allpermQ2` — permutation distribution of Q²  
- `results.Q2` — observed nested CV Q²  
- `results.allbootQ2` — bootstrap Q² distribution  
- `results.Q2_CI` — bootstrap confidence interval  
- `results.learningSizes` — sample sizes evaluated  
- `results.learningQ2` — Q² values across sample sizes  

---

## `XL` (numeric, [p×LV])

X-loadings matrix from `plsregress`, used to map LV loadings to ROIs.

If `XL` is empty, the function attempts to use:

```
results.finalXLoadings
```

**Important:** `XL` must correspond to the same feature ordering as the columns in `X`.

---

## `roiNames` (cellstr/string, length p)

ROI/feature labels.

If empty, defaults to:

```
ROI_001 ... ROI_p
```

---

## `atlasFile` (char/string)

Path to a labeled atlas NIfTI used for ROI mapping.

Atlas requirements:

- integer labels `1..p`
- `0` used for background

Mapping rule:

```
feature i → voxels where atlasData == i
```

---

# Name–Value Options

| Option | Default | Meaning |
|---|---:|---|
| `TopN` | 20 | Number of rows to display at top of ROI table (table file exports all rows). |
| `TopK` | 20 | Number of top features (by absolute meanFeatureWeight / meanBeta) for bar/table export. |
| `LV` | 1 | LV index to visualize as a brain map. |
| `OutPrefix` | `'PLSR'` | Prefix for exported files (CSV + NIfTI). |
| `VIP_thresh` | 1 | VIP threshold for “robust contributor” marking/labeling. |
| `stab_thresh` | 2 | stabilityZ threshold for “robust contributor” marking/labeling. |
| `MapPrctile` | 70 | Threshold percentile for LV map (keeps top `100−MapPrctile`% by \|loading\|). |
| `RelaxIfEmpty` | true | If no ROI passes robust thresholds, relax thresholds (75th percentile) for visualization only. |

---

# Outputs

## `ROI_table` (MATLAB table)

Returned table with columns:

- `ROI` — ROI/feature label  
- `VIP` — VIP score  
- `stabilityZ` — stabilityZ statistic  
- `RobustContributor` — boolean based on thresholds  

Optional columns if available:

- `meanBeta`  
- `selectionFrequency`  
- `signStability`  
- `featureStability`

The table is also saved to disk as:

```
'<OutPrefix>_ROI_VIP_stability.csv'
```

---

# Files created (exports)

## CSV

```
'<OutPrefix>_ROI_VIP_stability.csv'
```

Full ROI table (VIP, stabilityZ, robust flag).

If weights are available:

```
'<OutPrefix>_top<k>_weights.csv'
```

Top-K weights table.

---

## NIfTI maps (requires SPM)

```
'<OutPrefix>_VIP_map.nii'
```

VIP values mapped to atlas.

```
'<OutPrefix>_stabilityZ_map.nii'
```

stabilityZ mapped to atlas.

```
'<OutPrefix>_LV<LV>_map.nii'
```

LV loadings mapped to atlas (normalized to max \|loading\|).

```
'<OutPrefix>_LV<LV>_map_thresh.nii'
```

Thresholded LV map based on `MapPrctile`.

---

# Figures created (interactive MATLAB figures)

## 1. Scatter plot: VIP vs stabilityZ

- X-axis: VIP  
- Y-axis: stabilityZ  
- Horizontal/vertical lines at `VIP_thresh` and `stab_thresh`  
- Robust ROIs labeled with text.

---

## 2. Multi-slice 3×N panel

Rows:

1. LV map (selected LV)  
2. VIP map (thresholded at `VIP_thresh`)  
3. stabilityZ map (thresholded at `stab_thresh`)

Robust ROI labels are overlaid on each slice.

---

## 3. Top-K bar plot (optional)

Requires `results.meanFeatureWeight` or `results.meanBeta`.

Bars show the **mean CV weight** for the top-K absolute weights.

---

## 4. Feature stability stem plot (optional)

Requires:

```
results.featureStability
```

Shows proportion of runs where the feature weight is non-zero.

---

## 5. Predicted vs observed scatter (optional)

Requires:

```
results.cvObserved
results.cvPredicted
```

The plot includes:

- identity line  
- OLS regression line  
- robust regression line  
- optional density shading  
- optional coloring by subject across repeats  
- optional subject-trajectory lines across repeats  

This visualization helps identify:

- systematic prediction bias  
- unstable predictions across repeats  
- subject-specific prediction patterns

---

## 6. Permutation Q² histogram (optional)

Requires:

```
results.allpermQ2
```

Displays:

- permutation distribution  
- vertical line at observed Q²

Used to evaluate whether predictive performance exceeds chance.

---

## 7. Bootstrap Q² histogram (optional)

Requires:

```
results.allbootQ2
```

Displays the out-of-bag bootstrap distribution of Q².

---

## 8. Learning curve

Requires:

```
results.learningSizes
results.learningQ2
```

Plots predictive performance as a function of training sample size.

---

# Robust contributor definition

By default:

```
RobustContributor = (VIP > VIP_thresh) & (abs(stabilityZ) > stab_thresh)
```

This mirrors a common PLS reporting convention:

- **VIP** captures importance in the fitted model  
- **stabilityZ** captures robustness across cross-validated model refits  

---

## When `RelaxIfEmpty = true`

If **no** ROI passes thresholds, the function relaxes them **for visualization only**:

```
VIP_thresh_vis  = prctile(VIP,75)
stab_thresh_vis = prctile(abs(stabilityZ),75)
```

The exported table reflects these thresholds in that scenario.

To enforce strict thresholds:

```
RelaxIfEmpty = false
```

---

# Atlas / feature alignment checklist

This function assumes:

1. The **column order** of feature matrix `X` matches  
2. The **row order** of `results.VIP`, `results.stabilityZ`, and `XL`, and  
3. Atlas labels `1..p` correspond to those same features.

Common failure modes:

- atlas uses non-contiguous labels (e.g., 1001..2000)  
- atlas has fewer labels than `p`  
- ROI list includes background index

If your atlas labels are not `1..p`, remap them before using this function.

---

# Interpreting the LV map

LV loadings represent a **multivariate pattern** rather than independent univariate effects.

Key points:

- positive and negative loadings represent opposite sides of the latent pattern  
- LV maps should be interpreted as **pattern descriptors**, not isolated effects  
- sign interpretation may depend on the direction of the outcome variable

---

# Dependencies

Requires:

- **SPM** (for atlas I/O):  
  `spm_vol`, `spm_read_vols`, `spm_write_vol`

Optional:

- `robustfit` (Statistics and Machine Learning Toolbox) for robust regression line.

---

# Recommended use by study type

## ROI-level features (PET binding, fMRI betas, morphometry)

Atlas mapping is straightforward.

Typical thresholds:

```
VIP_thresh = 1
stab_thresh = 2
```

---

## Connectivity edges / graph metrics

If features represent **edges**, direct atlas mapping may not apply.

You can still use:

- ROI table  
- VIP vs stability scatter  
- top-K plots  

but skip atlas/NIfTI mapping.

---

## p/n ratio considerations

When **p is large and n small**:

- stabilityZ may be noisy  
- rely more on selectionFrequency and signStability diagnostics  

When **p is moderate and n larger**:

- stabilityZ becomes more interpretable  
- thresholds may be tightened.

---

# Minimal reproducible call

```matlab
ROI_table = plot_PLSR_diagnostics_neuroimaging(results, results.finalXLoadings, roiNames, atlasFile, ...
    'LV', 1, 'OutPrefix', 'Study1_PLSR', 'TopN', 20, 'TopK', 20);
```