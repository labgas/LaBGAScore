# TSPO_ENet_pipeline — User Guide

## Overview

`TSPO_ENet_pipeline` implements a robust **Elastic Net classification pipeline** for TSPO PET ROI feature matrices.

Elastic Net combines:

- **L1 (lasso) regularization** → feature selection
- **L2 (ridge) regularization** → shrinkage and stability

This allows the model to:

- handle **high-dimensional data**
- reduce **overfitting**
- provide **sparse and interpretable feature weights**
- remain useful when **p ≈ n or p > n**

The pipeline is designed primarily for **TSPO PET ROI datasets**, but the overall logic is also applicable to other neuroimaging feature matrices.

The architecture emphasizes:

- repeated nested cross-validation
- leakage-free preprocessing
- hyperparameter tuning
- permutation testing
- **out-of-bag (OOB) bootstrap confidence intervals**
- feature stability metrics
- learning curves

---

## 1. Pipeline Architecture

The pipeline uses **repeated nested cross-validation**.

### Outer Cross-Validation

Purpose: estimate **generalization performance**.

Workflow:

1. Split data into **K folds**
2. Train model on **K−1 folds**
3. Test on the **held-out fold**

Outputs:

- AUC
- accuracy
- sensitivity
- specificity

### Inner Cross-Validation

Purpose: tune Elastic Net hyperparameters:

- **alpha** – L1 vs L2 balance
- **lambda** – regularization strength

For each outer training fold:

1. Evaluate the alpha × lambda grid
2. Select the combination with the highest inner-fold AUC

### Repeated Cross-Validation

The outer CV procedure is repeated **nRepeats** times.

Benefits:

- reduces variance of performance estimates
- improves robustness of feature selection
- enables feature stability summaries

---

## 2. Input Data Structure

### Feature Matrix X

`X : [n × p]`

| Dimension | Meaning |
|---|---|
| n | subjects |
| p | features / ROIs |

Typical TSPO PET features:

- ROI binding estimates
- atlas-based regional uptake values
- other regional PET summary measures

### Outcome Vector Y

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

## 3. Preprocessing

Scaling is performed **inside each cross-validation fold**:

`Xtrain_z = (Xtrain − mean(Xtrain)) / std(Xtrain)`  
`Xtest_z  = (Xtest − mean(Xtrain)) / std(Xtrain)`

This prevents leakage from the held-out test data into training.

For TSPO ROI features, z-scoring is generally appropriate because ROIs can differ in scale and variance.

---

## 4. Hyperparameters

Default settings are tuned for **small-sample neuroimaging studies**.

| Parameter | Default | Purpose |
|---|---|---|
| outerK | 5 | outer CV folds |
| innerK | 4 | inner CV folds |
| nRepeats | 50 | repeated outer CV |
| nPerm | 1000 | permutations |
| nBoot | 500 | bootstrap samples |
| learningSteps | 6 | learning curve points |

Elastic Net hyperparameters:

| Parameter | Default | Meaning |
|---|---|---|
| alphaGrid | [0.05 0.1 0.25 0.5 0.75 0.9 1] | lasso–ridge mixing |
| lambdaGrid | logspace(-3,1,25) | regularization strength |

Interpretation:

- **alpha = 1** → lasso-like, sparse model
- **alpha ≈ 0** → ridge-like, denser model
- **larger lambda** → stronger shrinkage

---

## 5. Output Structure

The pipeline returns a MATLAB structure:

`results`

### Cross-validated performance

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

### Hyperparameter selection

- `results.selectedAlpha`
- `results.selectedLambda`

### Model coefficients

- `results.betaStore`
- `results.interceptStore`
- `results.featureWeights`
- `results.meanFeatureWeight`

Elastic Net produces **sparse solutions**, so many coefficients may be exactly zero.

### Feature stability metrics

- `results.featureStability`
- `results.signStability`
- `results.selectionFrequency`

Interpretation guideline:

| Stability | Interpretation |
|---|---|
| >0.8 | highly stable |
| 0.4–0.8 | moderate |
| <0.4 | unstable |

---

## 6. Global Baseline Model

The pipeline also fits a simple baseline model using:

`mean(X,2)`

This tests whether a global TSPO PET signal alone can classify the groups.

Output:

`results.AUC_global`

Interpretation example:

- `AUC_global = 0.58`
- `AUC_model = 0.72`

This suggests regional patterning adds predictive information beyond a global signal shift.

---

## 7. Interpreting Elastic Net Coefficients

Positive coefficient → higher feature values in the positive class  
Negative coefficient → lower feature values in the positive class

Example:

- insula weight = +0.25
- precuneus weight = −0.19

Interpretation: the positive group tends to show relatively higher insula binding and lower precuneus binding.

Because coefficients are estimated in a multivariate penalized model, interpretation should focus on:

- **stability across resampling**
- **selection frequency**
- **sign consistency**

rather than a single coefficient from a single fold.

---

## 8. Permutation Testing

Permutation testing evaluates whether predictive performance exceeds chance.

Procedure:

1. Shuffle labels
2. Re-run the fast cross-validated ENet routine
3. Compute AUC

Outputs:

- `results.allpermAUC`
- `results.permAUC`
- `results.permutation_p`

Interpretation:

| p-value | Meaning |
|---|---|
| <0.05 | significant |
| <0.01 | strong evidence |
| <0.001 | very strong evidence |

A healthy null distribution should be approximately symmetric and centered near **0.5**.

---

## 9. Bootstrap Confidence Intervals

Uncertainty in AUC is estimated using **out-of-bag (OOB) bootstrap sampling**.

Procedure:

1. Sample subjects **with replacement**
2. Train the model on the **in-bag sample**
3. Tune alpha/lambda **within the in-bag sample only**
4. Evaluate performance on the **out-of-bag subjects**
5. Repeat `nBoot` times

Outputs:

- `results.allbootAUC`
- `results.bootAUC`
- `results.AUC_CI`

Because evaluation is performed on **OOB subjects rather than the bootstrap sample itself**, this bootstrap estimate is usually **more conservative** than naive bootstrap performance estimates.

Accordingly:

- OOB bootstrap AUC may be **similar to or somewhat lower than** nested CV AUC
- nested CV remains the **primary estimate of generalization performance**
- bootstrap is best viewed as a **sampling variability summary**

---

## 10. Learning Curves

Learning curves estimate AUC as a function of sample size.

Outputs:

- `results.learningSizes`
- `results.learningAUC`

Interpretation:

- rising curve → more subjects would likely improve performance
- plateauing curve → model may be approaching its achievable limit

---

## 11. Typical Result Pattern

Example:

- Nested CV AUC = 0.70
- Permutation p = 0.02
- OOB bootstrap CI = [0.58–0.79]

Interpretation:

- classification is above chance
- permutation testing supports statistical significance
- bootstrap indicates uncertainty around performance
- bootstrap may be slightly lower than nested CV without implying an error

---

## 12. Recommended Reporting

Example manuscript text:

> Elastic Net classification with repeated nested cross-validation yielded AUC = 0.70.  
> Permutation testing (5000 permutations) confirmed performance exceeded chance (p = 0.02).  
> Out-of-bag bootstrap sampling indicated AUC confidence intervals of [0.58–0.79].  
> Regions with high selection frequency and stable coefficient signs were considered robust contributors.

---

## 13. Common Pitfalls

**Data leakage**  
Do not scale or select features on the full dataset before cross-validation.

**Class imbalance**  
If folds frequently miss one class, reduce `outerK` and/or `innerK`.

**Overinterpretation**  
Interpret stability metrics, not single-fold coefficients.

**Bootstrap interpretation**  
Bootstrap OOB AUC is not the primary performance estimate. Use nested CV AUC first, and bootstrap CI as a complement.

---

## 14. Summary

`TSPO_ENet_pipeline` provides:

- robust performance estimation
- automatic feature selection
- interpretable TSPO PET ROI weights
- permutation-based statistical inference
- out-of-bag bootstrap uncertainty estimation
- feature stability diagnostics

It is therefore well suited for **small-sample TSPO PET machine-learning studies**.
