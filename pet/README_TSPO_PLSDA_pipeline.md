# TSPO_PLSDA_pipeline — User Guide

## Overview

`TSPO_PLSDA_pipeline` implements a robust **Partial Least Squares Discriminant Analysis (PLS-DA)** pipeline for TSPO PET ROI feature matrices.

PLS-DA is particularly well suited to neuroimaging settings where:

- features are strongly correlated
- the number of features may approach or exceed the sample size
- prediction and interpretation both matter

The pipeline is designed for **subjects × features** datasets such as:

- TSPO PET ROI binding values
- atlas-based regional PET measures
- regional TSPO uptake estimates
- other TSPO ROI feature matrices

The architecture emphasizes:

- repeated nested cross-validation
- leakage-free preprocessing
- inner tuning of latent variables (LVs)
- permutation testing
- **out-of-bag (OOB) bootstrap confidence intervals**
- feature importance and stability metrics
- learning curves

---

# 1. Pipeline Architecture

The pipeline uses **repeated nested cross-validation**.

## Outer Cross-Validation

Purpose: estimate **generalization performance**.

Workflow:

1. Split data into **K folds**
2. Train on **K−1 folds**
3. Test on the **held-out fold**

Outputs include:

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

# 2. Input Data Structure

## Feature Matrix X

`X : [n × p]`

| Dimension | Meaning |
|---|---|
| n | subjects |
| p | ROI features |

Typical examples:

- TSPO PET ROI binding
- atlas-based PET regional summaries
- regional TSPO uptake measures

---

## Outcome Vector Y

`Y : [n × 1]`

Accepted formats:

- numeric
- logical
- categorical
- string
- cell array of strings

Internally converted to:

`yNum = double(Y == max(Y))`

The **maximum label becomes the positive class**.

---

# 3. Preprocessing

Scaling is performed **inside each cross-validation fold** to prevent leakage.

Example (z-score scaling):

`Xtrain_z = (Xtrain − mean(Xtrain)) / std(Xtrain)`  
`Xtest_z  = (Xtest − mean(Xtrain)) / std(Xtrain)`

For most TSPO ROI feature matrices, **z-score scaling is recommended**.

---

# 4. Hyperparameters

Default settings are designed for **small neuroimaging samples**.

| Parameter | Default | Purpose |
|---|---|---|
| outerK | 5 | outer CV folds |
| innerK | 4 | LV tuning |
| nRepeats | 50 | repeated outer CV |
| maxLV | 4 | maximum candidate latent variables |
| nPerm | 1000 | permutations |
| nBoot | 500 | bootstrap samples |
| learningSteps | 6 | learning curve points |

The actual LV cap is constrained by **sample size and matrix rank** within each training fold.

---

# 5. Output Structure

The pipeline returns a MATLAB structure:

`results`

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

Primary performance estimate:

`results.AUC`

---

## LV Selection and Model Weights

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

# 6. Global Baseline Model

A simple logistic model is also fitted using:

`mean(X,2)`

Output:

`results.AUC_global`

This helps determine whether a **distributed TSPO regional pattern** improves prediction beyond a global signal shift.

---

# 7. Interpreting PLS-DA Outputs

## Latent Variables

PLS-DA summarizes covariance between ROI features and outcome into **latent variables (LVs)**.

Typical interpretation:

- **LV1** often reflects the strongest global TSPO pattern
- later LVs may explain less X variance but still contribute to discrimination

---

## VIP Scores

VIP reflects variable importance in the final model.

Common heuristic:

| VIP | Interpretation |
|---|---|
| <0.8 | weak |
| 0.8–1.0 | moderate |
| >1.0 | important |

---

## Stability Metrics

- `stabilityZ` → coefficient stability across CV runs
- `signStability` → sign consistency
- `selectionFrequency` → frequency of appearing among top-K weights

These metrics should be interpreted **together**, rather than relying on a single coefficient.

---

# 8. Permutation Testing

Permutation testing evaluates whether predictive performance exceeds chance.

Procedure:

1. Shuffle labels
2. Re-run the quick cross-validated PLS-DA routine
3. Compute AUC

Outputs:

- `results.allpermAUC`
- `results.permAUC`
- `results.permutation_p`

A healthy null distribution should be **approximately centered near 0.5**.

---

# 9. Bootstrap Confidence Intervals

Uncertainty in AUC is estimated using **out-of-bag (OOB) bootstrap sampling**.

Procedure:

1. Sample subjects **with replacement**
2. Train the model on the **in-bag sample**
3. Tune latent variables within the in-bag sample
4. Evaluate AUC on the **out-of-bag subjects**
5. Repeat `nBoot` times

Outputs:

- `results.allbootAUC`
- `results.bootAUC`
- `results.AUC_CI`

Because evaluation occurs on **OOB subjects rather than the bootstrap sample**, bootstrap estimates are typically **more conservative** than naive bootstrap approaches.

Nested cross-validation remains the **primary performance estimate**.

---

# 10. Learning Curves

Learning curves estimate AUC as a function of sample size.

Outputs:

- `results.learningSizes`
- `results.learningAUC`

Interpretation:

- increasing curve → more data likely improves performance
- plateau → model approaching its achievable limit

---

# 11. Typical Result Pattern

Example:

- Nested CV AUC = 0.68
- Permutation p = 0.02
- OOB bootstrap CI = [0.57–0.76]

Interpretation:

- classification is above chance
- permutation testing supports statistical significance
- bootstrap reflects uncertainty around performance
- bootstrap may be slightly lower than nested CV without indicating a coding error

---

# 12. Recommended Reporting

Example manuscript text:

> Partial Least Squares Discriminant Analysis of TSPO PET ROI features with repeated nested cross-validation yielded AUC = 0.68.  
> Permutation testing (5000 permutations) confirmed performance exceeded chance (p = 0.02).  
> Out-of-bag bootstrap sampling indicated AUC confidence intervals of [0.57–0.76].  
> Features with high VIP, stable sign, and high cross-validated stability were considered robust contributors.

---

# 13. Common Pitfalls

**Data leakage**  
Do not scale or select features on the full dataset before cross-validation.

**Class imbalance**  
If folds frequently miss one class, reduce `outerK` and/or `innerK`.

**Too many latent variables**  
PLS components must be constrained by rank and sample size.

**Bootstrap interpretation**  
Bootstrap OOB AUC is **not the primary performance estimate**. Use nested CV AUC first, and bootstrap CI as a complement.

---

# 14. Summary

`TSPO_PLSDA_pipeline` provides:

- robust performance estimation
- latent-variable modeling suited for correlated TSPO PET ROI features
- permutation-based statistical inference
- out-of-bag bootstrap uncertainty estimation
- feature importance and stability diagnostics

This makes it well suited for **small-sample TSPO PET machine-learning studies** where **correlated predictors and p≈n or p>n** are common.