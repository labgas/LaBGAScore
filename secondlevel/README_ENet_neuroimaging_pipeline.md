
# ENet_neuroimaging_pipeline / TSPO_ENet_pipeline — User Guide

## Overview

`ENet_neuroimaging_pipeline` and `TSPO_ENet_pipeline` implement a robust **Elastic Net classification pipeline** for neuroimaging feature matrices.

Elastic Net combines **L1 (lasso)** and **L2 (ridge)** regularization, allowing the model to:

- perform **feature selection**
- handle **high-dimensional data**
- reduce overfitting when **p ≈ n or p > n**

The pipelines are designed for **subjects × features** datasets common in neuroimaging such as:

- PET ROI binding values
- fMRI ROI beta estimates
- cortical thickness or morphometry
- connectivity matrices or edge features
- graph metrics
- multimodal ROI feature sets

The architecture mirrors the PLSDA pipelines and estimates predictive performance while minimizing overfitting by combining:

- repeated nested cross-validation
- leakage-free preprocessing
- hyperparameter tuning
- permutation testing
- bootstrap confidence intervals
- feature stability metrics
- learning curves

The **TSPO_ENet_pipeline** is the original TSPO-specific version.  
The **ENet_neuroimaging_pipeline** is the **generic version** supporting different scaling modes and global signal definitions.

---

# 1. Pipeline Architecture

The pipeline uses **repeated nested cross-validation**.

## Outer Cross-Validation

Purpose: estimate **generalization performance**.

Process:

Data → split into K folds  
Train model on K−1 folds  
Test model on held-out fold

Outputs include:

- AUC
- accuracy
- sensitivity
- specificity

---

## Inner Cross-Validation

Purpose: select the optimal **Elastic Net hyperparameters**:

- **alpha** → L1 vs L2 regularization balance
- **lambda** → regularization strength

Process:

Training set → inner CV  
Evaluate alpha × lambda grid  
Select combination with highest mean AUC

---

## Repeated Cross-Validation

The outer CV procedure is repeated:

nRepeats times

Benefits:

- reduces variance of performance estimates
- enables feature stability estimation
- improves reliability of interpretation

---

# 2. Input Data Structure

## Feature Matrix X

X : [n × p]

| Dimension | Meaning |
|---|---|
| n | subjects |
| p | features |

Examples:

- PET ROI binding
- fMRI ROI beta estimates
- cortical thickness
- graph metrics
- connectivity strengths
- edge weights

---

## Outcome Vector Y

Binary outcome vector:

Y : [n × 1]

Accepted formats:

- numeric
- logical
- categorical
- string
- cellstr

Internally converted to:

yNum = double(Y == max(Y))

The **maximum label becomes the positive class**.

---

# 3. Preprocessing

The pipeline performs **leakage-free scaling** inside each cross-validation fold.

Example:

Xtrain_z = (Xtrain − mean(Xtrain)) / std(Xtrain)  
Xtest_z  = (Xtest − mean(Xtrain)) / std(Xtrain)

Scaling parameters are learned **only from training data**.

Scaling options:

- opts.scale = 'zscore' (default)
- opts.scale = 'center'
- opts.scale = 'none'

For most neuroimaging matrices **z-score scaling is recommended**.

---

# 4. Choosing Hyperparameters

Default parameters are optimized for **small neuroimaging samples**.

| Parameter | Default | Purpose |
|---|---|---|
outerK | 5 | outer CV folds |
innerK | 4 | hyperparameter tuning |
nRepeats | 50 | CV repetitions |
nPerm | 1000 | permutation tests |
nBoot | 500 | bootstrap samples |
learningSteps | 6 | learning curve points |

Elastic Net hyperparameters:

| Parameter | Default | Purpose |
|---|---|---|
alphaGrid | [0.05–1] | L1/L2 mixing |
lambdaGrid | logspace(-3,1,25) | regularization strength |

---

# 5. Recommended Settings by Dataset Type

## High-dimensional datasets (p >> n)

Examples:

- connectivity edges
- voxelwise features
- multimodal matrices

Recommended:

outerK = 4–5  
innerK = 3–4  
nRepeats ≥ 50

Higher alpha encourages **sparser models**.

---

## Moderate dimensional datasets (p ≈ n)

Examples:

- ROI beta estimates
- cortical parcellations

Recommended:

outerK = 5  
innerK = 4–5  
nRepeats = 30–50

---

## Lower dimensional datasets (p < n)

Examples:

- graph summary metrics
- small ROI feature sets

Recommended:

outerK = 5–10  
innerK = 5  
nRepeats = 10–30

---

# 6. Output Structure

The pipeline returns:

results

### Cross‑validated performance

results.AUC  
results.ACC  
results.SENS  
results.SPEC  

Fold‑level metrics:

results.allAUC  
results.allACC  
results.allSENS  
results.allSPEC  

Primary metric:

results.AUC

---

### Hyperparameter selection

results.selectedAlpha  
results.selectedLambda

---

### Model coefficients

results.betaStore  
results.featureWeights  
results.meanFeatureWeight

Elastic Net coefficients represent **feature importance**.

Due to L1 regularization many coefficients may be **exactly zero**.

---

### Feature stability metrics

results.featureStability  
results.signStability  
results.selectionFrequency

Interpretation guidelines:

| Stability | Meaning |
|---|---|
>0.8 | highly stable |
0.4–0.8 | moderate |
<0.4 | unstable |

---

# 7. Global Baseline Model

The pipeline also fits a simple baseline model using a global feature:

mean(X,2)

or

median(X,2)

depending on:

opts.globalFun

Example interpretation:

AUC_global = 0.60  
AUC_model = 0.74  

Regional patterns provide additional predictive information.

---

# 8. Interpreting Elastic Net Coefficients

Elastic Net produces **sparse linear models**.

Positive coefficient → higher values in positive class  
Negative coefficient → lower values in positive class

Example:

insula weight = +0.28  
precuneus weight = −0.21

Interpretation: patients show higher insula values and lower precuneus values.

---

# 9. Permutation Testing

Permutation testing evaluates whether performance exceeds chance.

Procedure:

shuffle labels → run CV → compute AUC

Output:

results.permutation_p

Interpretation:

| p-value | Meaning |
|---|---|
<0.05 | significant |
<0.01 | strong evidence |
<0.001 | very strong evidence |

---

# 10. Bootstrap Confidence Intervals

Bootstrap estimates uncertainty in AUC.

Procedure:

resample subjects with replacement → run quickCV → repeat nBoot times

Output:

results.AUC_CI

Example:

AUC = 0.72  
CI = [0.60 0.85]

Bootstrap AUC may be optimistic in small samples.  
Nested CV remains the **primary performance estimate**.

---

# 11. Learning Curves

Learning curves estimate **performance vs sample size**.

If curve increases:

more subjects likely improve performance.

If plateaued:

dataset near achievable performance.

---

# 12. Typical Result Pattern

Example:

Nested CV AUC = 0.74  
Permutation p = 0.01  
Bootstrap CI = [0.78–0.97]

Interpretation:

- model performs above chance
- bootstrap optimistic for small n
- nested CV is main estimate

---

# 13. Recommended Reporting

Example manuscript text:

Elastic Net classification with repeated nested cross-validation yielded AUC = 0.74.  
Permutation testing (5000 permutations) confirmed performance exceeded chance (p = 0.01).  
Bootstrap sampling indicated AUC 0.90 [0.78–0.97], reflecting optimism in small samples.  
Regions with high selection frequency and stable coefficient signs were considered robust contributors.

---

# 14. Common Pitfalls

### Data leakage

Never perform scaling or feature selection on the **full dataset before CV**.

### Class imbalance

If CV warnings occur reduce outerK (e.g., outerK = 3).

### Overinterpreting single models

Interpret **stability metrics**, not individual weights.

---

# 15. Relationship to Plotting Scripts

Plotting functions typically use:

- results.featureWeights
- results.meanFeatureWeight
- results.selectionFrequency
- results.signStability

Ensure feature order matches ROI labels or atlas indices.

---

# 16. Summary

The Elastic Net pipelines provide:

- robust performance estimation
- automatic feature selection
- interpretable neuroimaging features
- statistical validation
- feature stability diagnostics

making them well suited for **small‑sample neuroimaging machine‑learning studies**.
