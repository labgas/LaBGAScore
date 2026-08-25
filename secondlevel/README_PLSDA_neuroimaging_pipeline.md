# PLSDA_neuroimaging_pipeline — User Guide

## Overview

`PLSDA_neuroimaging_pipeline` implements a robust **Partial Least Squares Discriminant Analysis (PLS-DA)** pipeline for neuroimaging feature matrices.

PLS-DA is particularly well suited to neuroimaging settings where:

- features are strongly correlated  
- the number of features may approach or exceed the sample size  
- prediction and interpretation both matter  

The pipeline is designed for **subjects × features** datasets such as:

- PET ROI binding values  
- fMRI ROI beta estimates  
- cortical thickness or morphometry  
- connectivity or graph-derived features  
- multimodal ROI feature matrices  

The architecture emphasizes:

- repeated nested cross-validation  
- leakage-free preprocessing  
- inner tuning of latent variables (LVs)  
- permutation testing  
- **out-of-bag (OOB) bootstrap confidence intervals**  
- feature importance and stability metrics  
- learning curves  

---

# Pipeline Architecture

The pipeline uses **repeated nested cross-validation**.

## Outer Cross-Validation

Purpose: estimate **generalization performance**.

Workflow:

1. Split data into **K folds**
2. Train on **K−1 folds**
3. Test on the **held-out fold**

Metrics:

- AUC  
- accuracy  
- sensitivity  
- specificity  

---

## Inner Cross-Validation

Purpose: select the optimal number of **latent variables (LVs)**.

For each outer training fold:

1. Evaluate LV = 1 … maxLV  
2. Select the LV with the highest inner-CV AUC  

---

## Repeated Cross-Validation

The outer CV procedure is repeated **nRepeats** times.

Benefits:

- reduces variance of performance estimates  
- improves robustness of LV selection  
- enables feature stability summaries  

---

# Input Data Structure

## Feature Matrix

```
X : [n × p]
```

| Dimension | Meaning |
|-----------|--------|
| n | subjects |
| p | features |

Example features:

- PET ROI binding  
- fMRI ROI beta estimates  
- cortical thickness  
- graph metrics  
- connectivity strengths  

---

## Outcome Vector

```
Y : [n × 1]
```

Accepted formats:

- numeric  
- logical  
- categorical  
- string  
- cell array of strings  

Internally converted to:

```
yNum = double(Y == max(Y))
```

The **maximum label becomes the positive class**.

---

# Preprocessing

Scaling is performed **inside each cross-validation fold** to prevent leakage.

Example:

```
Xtrain_z = (Xtrain − mean(Xtrain)) / std(Xtrain)
Xtest_z  = (Xtest − mean(Xtrain)) / std(Xtrain)
```

Supported scaling modes:

```
opts.scale = 'zscore'   (default)
opts.scale = 'center'
opts.scale = 'none'
```

For most neuroimaging feature matrices **z-score scaling is recommended**.

---

## Covariate control (fold-wise nuisance regression)

Covariates are regressed out **inside every cross-validation fold**, with the
nuisance coefficients estimated on the training fold only and then applied to
both folds:

```
opts.covariates     = [n x nCov] numeric, one ROW per subject   (default [])
opts.covariateNames = {'age','sex'}                              (default auto)
```

Do **not** add a column of ones; the intercept is handled internally. Encode
categorical covariates as dummy columns beforehand. Leaving `opts.covariates`
empty reproduces the previous behaviour exactly.

The order of operations is **residualize, then scale**, and it is applied at
every preprocessing site: outer folds, inner tuning folds, bootstrap resamples
and learning-curve subsamples. Two consequences worth knowing:

- The scale constants are computed on the *residualized* data, so the model's
  inputs really do have unit variance in the space it operates in.
- Residualization lowers `rank(Xtrain)` by up to `rank(C)`, so the latent-
  variable cap is applied after it.

### Do not pre-residualize

Both calling scripts previously instructed users to residualize the feature
matrix before running the pipeline. That fits the nuisance model on all
subjects, so each test fold influences the coefficients later used to transform
it.

**The direction of the resulting bias is not the intuitive one.** Measured on
synthetic data where a covariate drives *all* the apparent association (so the
correct answer is chance):

| route | PLS-DA AUC (chance 0.5) | PLSR Q² (chance ≈ 0) |
|---|---|---|
| raw X, no control | 0.914 | +0.730 |
| pre-residualized full matrix | **0.332** | **−0.766** |
| fold-wise (`opts.covariates`) | **0.525** | **−0.081** |

Full-sample residualization lands far *below* chance, not above: over-removal
subtracts a component estimated partly from the test fold, which inverts the
test-fold relationship. This is the same effect reported for confound
regression in decoding analyses. Where X carries genuine signal beyond the
confound the two routes agree closely (0.720 vs 0.731).

So the case for the fold-wise version is that it is **unbiased and
leakage-free**, returning the right answer in both regimes — not that it
removes optimism.

### Reproducibility

```
opts.seed = 1    (default)
```

Results are now reproducible from run to run **and independent of parallel pool
size**. Previously the permutation, bootstrap and learning-curve stages drew
their randomness on workers, whose streams the client seed never reached.

---

### Note for classification

Y stays binary, so only X is residualized (`opts.residualizeY` is not offered —
a residualized Y would break `perfcurve`, `lassoglm` and the `scores > 0.5`
rule).

If a covariate is itself associated with the class, removing it from X also
removes part of the class signal, so the controlled estimate is **conservative
by design**. That is a property of confound control, not a defect.

---

# Hyperparameters

Default settings are designed for **small neuroimaging samples**.

| Parameter | Default |
|----------|---------|
| outerK | 5 |
| innerK | 4 |
| nRepeats | 50 |
| maxLV | 4 |
| nPerm | 1000 |
| nBoot | 500 |
| learningSteps | 6 |

The actual LV cap is constrained by **sample size and matrix rank** within each training fold.

---

# Output Structure

The pipeline returns a MATLAB structure:

```
results
```

## Cross-validated performance

- `results.AUC`
- `results.AUC_PR` (precision-recall AUC; useful when positives are rare)
- `results.ACC`
- `results.SENS`
- `results.SPEC`
- `results.ACC_balanced` (mean of sensitivity and specificity; useful under imbalance)

Fold-level metrics:

- `results.allAUC`
- `results.allAUC_PR`
- `results.allACC`
- `results.allSENS`
- `results.allSPEC`
- `results.allACC_balanced`

Primary estimate:

```
results.AUC
```

---

## Model Selection

- `results.selectedLV`
- `results.betaStore`
- `results.featureWeights`
- `results.meanFeatureWeight`

---

## Final Model (Interpretation Only)

- `results.finalLV`
- `results.betaFinal`
- `results.varExplainedX`
- `results.varExplainedY`
- `results.finalXLoadings`
- `results.finalYLoadings`

---

## Feature Importance / Stability

- `results.VIP`
- `results.meanBeta`
- `results.sdBeta`
- `results.stabilityZ`
- `results.signStability`

These metrics help identify **robust contributors across resampling**.

---

# Global Baseline Model

A simple logistic model using:

```
mean(X,2)
```

Output:

```
results.AUC_global
results.AUC_PR_global
```

This tests whether regional patterns outperform a global signal shift.

The baseline above is computed **in-sample**: the model is fit on all subjects
and predicts those same subjects, so it is not comparable to the
cross-validated model performance sitting next to it. A cross-validated
counterpart is now reported alongside it, run through the same repeated outer
CV:

```
results.AUC_global_cv
results.AUC_PR_global_cv
```

The in-sample fields are unchanged and retained for continuity. **Compare the
`_cv` fields against the model**, not the in-sample ones.

---


---

# Permutation Testing

Permutation testing evaluates whether predictive performance exceeds chance.

Procedure:

1. Shuffle labels  
2. Re-run the quick cross-validated PLS-DA routine  
3. Compute AUC (and precision-recall AUC)  

Outputs:

```
results.allpermAUC
results.permAUC
results.permutation_p
results.allpermAUC_PR
results.permAUC_PR
results.permutation_p_PR
```

A healthy null distribution should be **centered near 0.5**.

---

## The p-value uses a matched estimator

`results.permutation_p` is computed against `results.quickCV_observed`, **not**
against the headline `results.AUC`.

The reason: the null distribution is produced by the quick CV routine, while the
headline number comes from repeated nested CV. Those are two different
estimators, and the quick routine previously also ran untuned (maximum latent
variables), which overfits permuted data much harder than the tuned estimator
overfits real data. The null therefore sat far below the observed value.

Measured on true-null data, that mismatch produced a **false-positive rate near
40 % at alpha = 0.05** (mean p 0.10 where it should be about 0.5). Matching the
estimator on both sides restores calibration; the quick routine now also tunes
internally, which additionally recovers power.

`results.AUC` remains the number to report as performance. The p-value uses
the `(sum(perm >= obs) + 1) / (nValid + 1)` convention, so it can never be
exactly zero.

**p-values from runs predating this change are inflated.**

---

# Bootstrap Confidence Intervals

Uncertainty in AUC is estimated using **out-of-bag (OOB) bootstrap sampling**.

Procedure:

1. Sample subjects **with replacement**  
2. Train model on the **in-bag sample**  
3. Tune latent variables within the in-bag sample  
4. Evaluate AUC on the **out-of-bag subjects**  
5. Repeat `nBoot` times  

Outputs:

```
results.allbootAUC
results.bootAUC
results.AUC_CI
```

Because evaluation occurs on **OOB subjects rather than the bootstrap sample**, bootstrap estimates are typically **more conservative** than naive bootstrap approaches.

Nested cross-validation remains the **primary performance estimate**.

---

# Learning Curves

Outputs:

```
results.learningSizes
results.learningAUC
```

Interpretation:

- increasing curve → more data likely improves performance  
- plateau → model approaching its achievable limit  

---

# Summary

`PLSDA_neuroimaging_pipeline` provides:

- robust performance estimation  
- latent-variable modeling suited for correlated neuroimaging features  
- permutation-based statistical inference  
- OOB bootstrap uncertainty estimation  
- feature importance and stability diagnostics  

This makes it well suited for **small-sample neuroimaging machine-learning studies** where **p ≈ n or p > n** are common.