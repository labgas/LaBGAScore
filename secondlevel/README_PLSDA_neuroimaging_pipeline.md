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
- `results.ACC`
- `results.SENS`
- `results.SPEC`

Fold-level metrics:

- `results.allAUC`
- `results.allACC`
- `results.allSENS`
- `results.allSPEC`

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
- `results.selectionFrequency`

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
```

This tests whether regional patterns outperform a global signal shift.

---

# Permutation Testing

Permutation testing evaluates whether predictive performance exceeds chance.

Procedure:

1. Shuffle labels  
2. Re-run the quick cross-validated PLS-DA routine  
3. Compute AUC  

Outputs:

```
results.allpermAUC
results.permAUC
results.permutation_p
```

A healthy null distribution should be **centered near 0.5**.

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