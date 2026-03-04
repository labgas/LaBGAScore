# plot_ENet_diagnostics_neuroimaging --- User Guide

`plot_ENet_diagnostics_neuroimaging` is a visualization and diagnostics
helper for **Elastic Net neuroimaging pipelines**.

It mirrors the structure of the PLS‑DA plotting utilities and provides
standardized figures, tables, and brain maps for interpreting Elastic
Net models trained on ROI‑level neuroimaging data (PET, fMRI, structural
MRI, connectivity, etc.).

The function generates:

-   ROI importance tables
-   importance vs robustness scatter plots
-   NIfTI maps for model weights and stability
-   multi‑slice brain figures
-   top‑K feature plots
-   feature stability plots
-   hyperparameter selection histograms
-   optional post‑selection logistic regression statistics

The function assumes the atlas labels **1..p correspond to the feature
order in X**.

------------------------------------------------------------------------

# Quick start

``` matlab
% After running the ENet pipeline:
% results = ENet_neuroimaging_pipeline(X,Y,opts);

ROI_table = plot_ENet_diagnostics_neuroimaging(...
    results, X, Y, roiNames, atlasFile, ...
    'OutPrefix','ENet_run1');
```

------------------------------------------------------------------------

# Inputs

## results (struct)

Elastic Net results structure produced by the pipeline.

Required fields

    results.meanFeatureWeight   [p×1]
    results.featureWeights      [p×runs]

Optional fields

    results.selectionFrequency
    results.featureStability
    results.selectedAlpha
    results.selectedLambda

Descriptions

  -----------------------------------------------------------------------
  field                               meaning
  ----------------------------------- -----------------------------------
  meanFeatureWeight                   mean beta across CV runs

  featureWeights                      beta estimates across runs

  selectionFrequency                  proportion of runs where feature
                                      appears in Top‑K

  featureStability                    proportion of runs where beta ≠ 0

  selectedAlpha                       alpha selected in nested CV

  selectedLambda                      lambda selected in nested CV
  -----------------------------------------------------------------------

------------------------------------------------------------------------

## X

Feature matrix `[n × p]` used in the ENet pipeline.

Used for:

-   dimension checks
-   optional post‑selection logistic regression

------------------------------------------------------------------------

## Y

Binary outcome vector `[n × 1]`.

Accepted formats:

-   numeric
-   logical
-   categorical
-   string
-   cellstr

Internally converted to

    yNum = double(Y == max(Y))

------------------------------------------------------------------------

## roiNames

Vector of ROI labels.

If empty the function generates

    ROI_001
    ROI_002
    ...
    ROI_p

------------------------------------------------------------------------

## atlasFile

Path to labeled atlas NIfTI.

Requirements

    atlas labels 1..p = ROIs
    0 = background

Each voxel with label i receives the value of feature i.

------------------------------------------------------------------------

# Options

  option            default   description
  ----------------- --------- ----------------------------------------------
  TopN              20        rows printed from ROI table
  TopK              20        number of features in plots
  FreqThresh        0.5       robustness threshold for selection frequency
  WeightThresh      0         threshold for
  MapPrctile        70        percentile threshold for weight map
  DoPostSelection   true      perform logistic refit
  OutPrefix         ENet      prefix for exported files
  RelaxIfEmpty      true      relax thresholds if none pass

------------------------------------------------------------------------

# Output

## ROI_table

Returned MATLAB table

    ROI
    meanWeight
    absMeanWeight
    selectionFrequency
    featureStability
    RobustContributor

Also written to

    <OutPrefix>_ROI_weights_stability.csv

------------------------------------------------------------------------

# Files written

Tables

    <OutPrefix>_ROI_weights_stability.csv
    <OutPrefix>_postselection_refit_table.csv

NIfTI maps

    <OutPrefix>_meanWeight_map.nii
    <OutPrefix>_selectionFreq_map.nii
    <OutPrefix>_featureStability_map.nii
    <OutPrefix>_meanWeight_map_thresh.nii

------------------------------------------------------------------------

# Figures produced

## 1. Importance vs robustness

Scatter plot of

    |meanWeight| vs selectionFrequency

Features passing robustness thresholds are labeled.

------------------------------------------------------------------------

## 2. Brain slice visualization

Three rows:

    meanWeight map
    selectionFrequency map
    featureStability map

Slices are selected automatically from atlas regions.

------------------------------------------------------------------------

## 3. Top‑K feature weights

Bar plot of the strongest mean ENet weights.

------------------------------------------------------------------------

## 4. Feature stability plot

Stem plot of

    featureStability

representing the proportion of CV runs where the coefficient was
non‑zero.

------------------------------------------------------------------------

## 5. Hyperparameter histograms

If available

    results.selectedAlpha
    results.selectedLambda

histograms summarize the hyperparameters selected during nested CV.

------------------------------------------------------------------------

# Robust feature definition

Default rule

    |meanWeight| > WeightThresh
    AND
    selectionFrequency ≥ FreqThresh

Default values

    WeightThresh = 0
    FreqThresh = 0.5

If no feature passes thresholds and `RelaxIfEmpty=true`, thresholds are
relaxed to the **75th percentile** for visualization.

------------------------------------------------------------------------

# Optional post‑selection inference

If enabled:

1.  Top‑K non‑zero features are selected
2.  Logistic regression is fit

```{=html}
<!-- -->
```
    fitglm(Xsel, yNum, 'Distribution','binomial', ...
           'Link','logit','LikelihoodPenalty','jeffreys-prior')

Output table includes

    Beta_refit
    SE
    OddsRatio
    pValue

These statistics are **post‑selection summaries**, not formal inference.

------------------------------------------------------------------------

# Atlas alignment checklist

Ensure:

1.  columns of X correspond to features 1..p
2.  results.meanFeatureWeight uses the same order
3.  atlas labels correspond to those indices

Common issues

-   atlas labels not contiguous
-   atlas contains fewer labels than features
-   ROI order mismatch

------------------------------------------------------------------------

# Interpreting Elastic Net importance metrics

### meanFeatureWeight

Average coefficient across CV runs.

-   sign → direction of association
-   magnitude → relative importance

------------------------------------------------------------------------

### selectionFrequency

Proportion of runs where feature appears among top weights.

Typical interpretation

  frequency   meaning
  ----------- ---------------
  \>0.7       highly robust
  0.4--0.7    moderate
  \<0.4       unstable

------------------------------------------------------------------------

### featureStability

Proportion of runs where coefficient ≠ 0.

Reflects **sparsity stability** of the model.

------------------------------------------------------------------------

# Dependencies

Required:

-   MATLAB Statistics and Machine Learning Toolbox
-   SPM (spm_vol, spm_read_vols, spm_write_vol)

------------------------------------------------------------------------

# Minimal reproducible call

``` matlab
ROI_table = plot_ENet_diagnostics_neuroimaging(...
    results, X, Y, roiNames, atlasFile, ...
    'OutPrefix','Study1_ENet', ...
    'TopK',20, ...
    'FreqThresh',0.5);
```
