
# PLSDA_neuroimaging_pipeline — User Guide

## Overview

`PLSDA_neuroimaging_pipeline` implements a robust **Partial Least Squares Discriminant Analysis (PLS-DA)** classification pipeline for neuroimaging feature matrices.

It is designed for **subjects × features** datasets common in neuroimaging such as:

- PET ROI binding values
- fMRI ROI beta estimates
- cortical thickness or morphometry
- connectivity matrices or edge features
- graph metrics
- multimodal ROI feature sets

The pipeline estimates predictive performance while minimizing overfitting by combining:

- repeated nested cross-validation
- leakage-free preprocessing
- hyperparameter tuning (latent variables)
- permutation testing
- bootstrap confidence intervals
- feature stability metrics
- learning curves

The architecture is intentionally kept **very close to the TSPO-specific pipeline**, while remaining fully generic.

This approach is particularly useful when:

p ≈ n   or   p > n

which is common in neuroimaging datasets.

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

Purpose: select the optimal **number of latent variables (LV)**.

Process:

Training set → inner CV  
Evaluate LV = 1..maxLV  
Select LV with highest mean AUC  

This prevents hyperparameter overfitting.

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

Meaning:

The **maximum label becomes the positive class**.

---

# 3. Preprocessing

The pipeline performs **leakage-free scaling** inside each cross-validation fold.

Example:

Xtrain_z = (Xtrain − mean(Xtrain)) / std(Xtrain)  
Xtest_z  = (Xtest − mean(Xtrain)) / std(Xtrain)

Scaling parameters are learned **only from training data**.

Optional scaling modes:

opts.scale = 'zscore'   (default)  
opts.scale = 'center'  
opts.scale = 'none'

For most neuroimaging feature matrices, **z-scoring is recommended**.

---

# 4. Choosing Hyperparameters

Default parameters are optimized for **small neuroimaging samples**.

| Parameter | Default | Purpose |
|---|---|---|
outerK | 5 | outer CV folds |
innerK | 4 | LV tuning |
nRepeats | 50 | CV repetitions |
maxLV | 4 | max latent variables |
nPerm | 1000 | permutation tests |
nBoot | 500 | bootstrap samples |
learningSteps | 6 | learning curve points |

---

# 5. Recommended Settings by Dataset Type

## High-dimensional datasets (p >> n)

Example:

- connectivity edges
- multimodal features

Recommended:

outerK = 4–5  
innerK = 3–4  
maxLV = 2–4  
nRepeats ≥ 50  

---

## Moderate dimensional datasets (p ≈ n)

Example:

- ROI betas
- cortical parcellations

Recommended:

outerK = 5  
innerK = 4–5  
maxLV = 4–8  
nRepeats = 30–50  

---

## Lower dimensional datasets (p < n)

Example:

- small ROI feature sets

Recommended:

outerK = 5–10  
innerK = 5  
maxLV = 6–12  
nRepeats = 10–30  

---

# 6. Output Structure

The pipeline returns a structure:

results

## Cross-validated performance

results.AUC  
results.ACC  
results.SENS  
results.SPEC  

Fold-level metrics:

results.allAUC  
results.allACC  
results.allSENS  
results.allSPEC  

---

## Model selection

results.selectedLV

LV selected for each outer fold.

---

## Model weights

results.betaStore  
results.featureWeights  
results.meanFeatureWeight  
results.featureStability  

---

## Baseline global model

results.AUC_global

Logistic regression using a global summary feature (mean or median of X).

Purpose: test whether classification is driven by **global signal differences**.

---

## Final model (interpretation only)

results.finalLV  
results.betaFinal  
results.finalXLoadings  
results.finalYLoadings  
results.varExplainedX  
results.varExplainedY  

Important:

This model is **not used for performance estimation**.

---

# 7. Feature Importance Metrics

The pipeline provides multiple complementary measures.

## VIP scores

results.VIP

Interpretation:

| VIP | Interpretation |
|---|---|
<0.8 | weak |
0.8–1 | moderate |
>1 | important |
>1.5 | strong contributor |

---

## Weight stability

results.meanBeta  
results.sdBeta  
results.stabilityZ  

Where:

stabilityZ = meanBeta / sdBeta

Interpretation:

| stabilityZ | Meaning |
|---|---|
<1 | unstable |
1–2 | moderate |
>2 | stable |
>3 | very stable |

Common robust feature criterion:

VIP > 1  
AND  
|stabilityZ| > 2  

---

## Sign stability

results.signStability

Measures how consistently the **direction of the effect** appears across CV runs.

Example:

0.95 → highly consistent sign  
0.5  → unstable sign  

---

## Selection frequency

results.selectionFrequency

Proportion of runs where a feature appears among the **Top-20 absolute weights**.

Typical interpretation:

| Frequency | Meaning |
|---|---|
>0.7 | highly robust |
0.4–0.7 | moderate |
<0.4 | unstable |

---

# 8. Statistical Validation

## Permutation testing

Tests whether classification exceeds chance.

Procedure:

shuffle labels  
repeat CV  
compute AUC  

Output:

results.permutation_p

Interpretation:

| p | Meaning |
|---|---|
<0.05 | significant |
<0.01 | strong evidence |
<0.001 | very strong evidence |

---

## Bootstrap confidence intervals

Bootstrap estimates uncertainty in AUC.

results.AUC_CI

Example:

AUC = 0.72  
CI = [0.58 0.85]

Important:

Bootstrap AUC may be **higher than nested CV AUC** because bootstrap samples reuse observations.

Nested CV remains the **primary performance estimate**.

---

# 9. Learning Curves

Learning curves estimate:

performance vs sample size

Outputs:

results.learningSizes  
results.learningAUC  

Interpretation:

- increasing curve → more subjects would likely improve performance
- plateau → model may have reached dataset limits

---

# 10. Interpreting Latent Variable Brain Maps

PLS components represent **multivariate patterns**, not isolated regions.

Key principles:

### 1. Patterns rather than single regions

PLS components capture **distributed networks**.

### 2. Relative interpretation

Loadings should be interpreted **relative to other regions**.

Example pattern:

insula ↑  
ACC ↑  
precuneus ↓  

Interpretation:

patients show increased insula/ACC and decreased precuneus.

---

### 3. Later LVs can still be important

It is common that:

- LV1 explains large variance in X
- LV2+ explain smaller variance in X but **more variance in Y**

These later components may capture **disease-specific patterns**.

---

# 11. Common Pitfalls

### Data leakage

Never perform:

- scaling
- feature selection
- nuisance regression involving Y

on the full dataset before cross-validation.

---

### Missing classes in folds

If warnings occur reduce:

outerK  
innerK  

---

### Overinterpreting single weights

Always interpret weights together with:

- VIP
- stabilityZ
- selectionFrequency

---

# 12. Summary

`PLSDA_neuroimaging_pipeline` provides:

- robust performance estimation
- interpretable multivariate brain patterns
- statistical validation
- feature stability diagnostics

making it well suited for **small-sample neuroimaging machine learning studies**.
