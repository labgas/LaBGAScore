# ENet Neuroimaging Pipeline

Robust **Elastic Net classification pipeline** for neuroimaging feature
matrices.

This repository provides two closely related MATLAB pipelines:

-   **ENet_neuroimaging_pipeline** -- generic implementation for
    neuroimaging feature matrices
-   **TSPO_ENet_pipeline** -- original TSPO PET--specific version

Both implement **Elastic Net regularized classification** using repeated
nested cross‑validation, permutation testing, and out‑of‑bag bootstrap
uncertainty estimation.

The pipelines are designed for **small to moderate sample size
neuroimaging datasets** where feature dimensionality can approach or
exceed sample size.

------------------------------------------------------------------------

# Overview

Elastic Net combines:

-   **L1 (lasso) regularization** → feature selection
-   **L2 (ridge) regularization** → shrinkage and stability

This combination allows models to:

-   handle **high‑dimensional data (p ≥ n)**
-   perform **automatic feature selection**
-   reduce **overfitting**
-   provide **interpretable feature weights**

Typical neuroimaging inputs include:

-   PET ROI binding values
-   fMRI ROI beta estimates
-   cortical thickness / morphometry
-   functional connectivity edges
-   graph metrics
-   multimodal ROI feature matrices

------------------------------------------------------------------------

# Pipeline Architecture

The pipelines estimate predictive performance using **repeated nested
cross‑validation**.

## Outer Cross‑Validation

Purpose: estimate **generalization performance**.

Workflow:

1.  Split data into **K folds**
2.  Train model on **K−1 folds**
3.  Test on the **held‑out fold**

Metrics:

-   AUC (primary metric)
-   accuracy
-   sensitivity
-   specificity

------------------------------------------------------------------------

## Inner Cross‑Validation

Purpose: tune Elastic Net hyperparameters:

-   **alpha** -- L1 vs L2 regularization balance
-   **lambda** -- regularization strength

For each outer training fold:

1.  Evaluate the **alpha × lambda grid**
2.  Select combination with **highest inner‑CV AUC**

------------------------------------------------------------------------

## Repeated Cross‑Validation

The outer CV procedure is repeated **nRepeats times**.

Benefits:

-   reduces variance of performance estimates
-   enables feature stability analysis
-   improves reliability of interpretation

------------------------------------------------------------------------

# Input Data Structure

## Feature Matrix

X : \[n × p\]

  Dimension   Meaning
  ----------- --------------------
  n           number of subjects
  p           number of features

Example features:

-   PET ROI binding
-   fMRI ROI beta estimates
-   cortical thickness
-   connectivity strengths
-   graph metrics

------------------------------------------------------------------------

## Outcome Vector

Y : \[n × 1\]

Accepted formats:

-   numeric
-   logical
-   categorical
-   string
-   cell array of strings

Internally converted to:

yNum = double(Y == max(Y))

The **maximum label becomes the positive class**.

------------------------------------------------------------------------

# Preprocessing

Scaling is performed **inside each cross‑validation fold** to prevent
data leakage.

Example (z‑score):

Xtrain_z = (Xtrain − mean(Xtrain)) / std(Xtrain)\
Xtest_z = (Xtest − mean(Xtrain)) / std(Xtrain)

Scaling options:

opts.scale = 'zscore' (default)\
opts.scale = 'center'\
opts.scale = 'none'

For most neuroimaging datasets **z‑score scaling is recommended**.

------------------------------------------------------------------------

# Hyperparameters

Default settings are designed for **small neuroimaging datasets**.

  Parameter       Default
  --------------- ---------
  outerK          5
  innerK          4
  nRepeats        50
  nPerm           1000
  nBoot           500
  learningSteps   6

Elastic Net parameters:

  Parameter    Default
  ------------ ----------------------------------
  alphaGrid    \[0.05 0.1 0.25 0.5 0.75 0.9 1\]
  lambdaGrid   logspace(-3,1,25)

------------------------------------------------------------------------

# Output Structure

The pipeline returns a MATLAB structure:

results

## Cross‑validated performance

-   results.AUC
-   results.ACC
-   results.SENS
-   results.SPEC

Fold‑level metrics:

-   results.allAUC
-   results.allACC
-   results.allSENS
-   results.allSPEC

The **primary performance estimate** is:

results.AUC

------------------------------------------------------------------------

## Hyperparameter Selection

results.selectedAlpha\
results.selectedLambda

------------------------------------------------------------------------

## Model Coefficients

results.betaStore\
results.featureWeights\
results.meanFeatureWeight

Elastic Net produces **sparse models**, meaning many coefficients may be
exactly zero.

------------------------------------------------------------------------

## Feature Stability Metrics

results.featureStability\
results.signStability\
results.selectionFrequency

Interpretation guideline:

  Stability   Interpretation
  ----------- ----------------
  \>0.8       highly stable
  0.4--0.8    moderate
  \<0.4       unstable

------------------------------------------------------------------------

# Global Baseline Model

A simple baseline logistic model is also computed using a **global
feature**:

mean(X,2)

or

median(X,2)

This helps determine whether **regional information improves
prediction** beyond global signal.

------------------------------------------------------------------------

# Permutation Testing

Permutation testing evaluates whether performance exceeds chance.

Procedure:

1.  Shuffle outcome labels
2.  Re‑run cross‑validation
3.  Compute AUC

Output:

results.permutation_p

------------------------------------------------------------------------

# Bootstrap Confidence Intervals

Uncertainty in AUC is estimated using **out‑of‑bag (OOB) bootstrap
sampling**.

Procedure:

1.  Sample subjects **with replacement**
2.  Train model on the **in‑bag sample**
3.  Evaluate performance on **out‑of‑bag subjects**
4.  Repeat nBoot times

Outputs:

results.allbootAUC\
results.bootAUC\
results.AUC_CI

Because evaluation occurs on **OOB subjects rather than the bootstrap
sample**, bootstrap estimates are typically **more conservative** than
naive bootstrap approaches.

Nested cross‑validation remains the **primary performance estimate**.

------------------------------------------------------------------------

# Learning Curves

Learning curves estimate how performance changes with **sample size**.

results.learningSizes\
results.learningAUC

------------------------------------------------------------------------

# Summary

The Elastic Net pipelines provide:

-   robust performance estimation
-   automatic feature selection
-   interpretable neuroimaging features
-   permutation‑based statistical inference
-   bootstrap uncertainty estimation
-   feature stability diagnostics

These properties make the pipelines well suited for **small‑sample
neuroimaging machine‑learning studies**.
