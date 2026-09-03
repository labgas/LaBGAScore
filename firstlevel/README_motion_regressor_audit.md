# Motion regressor audit — state of play

Working note, 2026/09/03. Written as a handover: enough to pick this up cold, on
another machine, without the conversation it came from.

The short version: the first-level motion regressors have never been what the
script headers say they are. Whether that matters turns out to depend entirely
on the study, and of the three checked so far, one needs a re-analysis, one does
not, and one is clean.

## Contents

- [What was found](#what-was-found)
- [How much it matters, in principle](#how-much-it-matters-in-principle)
- [The two tools](#the-two-tools)
- [Results per study](#results-per-study)
- [Verdicts](#verdicts)
- [What is settled and what is not](#what-is-settled-and-what-is-not)
- [Repository state](#repository-state)
- [Next steps](#next-steps)

## What was found

While writing [`README.md`](README.md) and [`README_phMRI.md`](README_phMRI.md),
reading the scripts end to end surfaced a set of defects. Most were fixed on
`main`; one is held back on the `firstlevel-fixes` branch pending this audit.

In `s2_fit_model`, `s2a_fit_model_multisess_multitask` (twice) and
`s1b_fit_phMRI_model`, the motion block reads:

```matlab
zscore = @(x) zscore(x);                      % defined
Rmotionzscore = varfun(deriv, Rmotionderiv);  % but never used - deriv again
```

The wrapper is defined and then not called, so the step that should z-score
applies `gradient` a second time. The model therefore contains the **first and
second derivatives** of the realignment parameters, where the headers describe
"six directions and their derivatives, z-scored", and where the CANlab/Friston
24-parameter model expects the parameters themselves plus their derivatives.

Measured on a real confound file: the first output column is bit-identical to
`gradient(trans_x)` and correlates **−0.010** with `trans_x`; under the fix it
has mean 0, SD 1 and correlates **1.000** with it. The two 12-column blocks
share about **64%** of their variance, so the model is not missing the
information entirely, but it is regressing out a substantially different
subspace. The sustained-displacement directions are the ones absent, and the
180 s high-pass filter does not recover them (subspace overlap 0.54 raw, 0.58
after filtering).

**This is not documented anywhere in the repo, deliberately.** It was kept out
of the READMEs and script comments pending internal discussion.

## How much it matters, in principle

Simulations on synthetic data with a known answer, before any real data was
touched. The estimator is not biased by the wrong parameterisation as such:

| Regime | Bias in β | SD of β | Group false-positive rate |
|---|---|---|---|
| Motion independent of task | ~0 | 0.058 → 0.167 | unaffected |
| Task-locked, artifact direction random across subjects | ~0 | 0.066 → 0.128 | unaffected |
| Task-locked, artifact **consistent** across subjects | 0 → +0.024 | — | 4.3% → **23.7%** |

The first two regimes cost power and nothing else: the extra scatter lands in
the between-subject variance and group inference stays valid. Only the third
biases a group map, and it requires the motion to be both task-locked **and**
pushed the same way in most subjects.

That distinction is what the diagnostic below measures, and it is the thing to
look at first for any study.

## The two tools

Both in `firstlevel/functions/`, both on `main`.

**`LaBGAScore_firstlevel_task_motion_diagnostics`** — reads only the fMRIPrep
confound files and the BIDS events files, so it runs in seconds and needs no
`SPM.mat` and no betas. Reports, per condition, the correlation between the
convolved regressor and framewise displacement plus a one-sample *t* across
subjects (the consistency test), and per contrast a leakage bound: how far a
motion artifact can displace that contrast, in units of its own standard error,
per unit ratio of artifact to residual noise.

```matlab
[T,G] = LaBGAScore_firstlevel_task_motion_diagnostics(BIDSdir, fullfile(derivdir,'fmriprep'), ...
    'dsgn', DSGN, 'tr', DSGN.tr, 'hpf', DSGN.hpf, ...
    'quadratic', LaBGAS_options.movement_reg_quadratic, 'task', '<label>')
```

**`LaBGAScore_firstlevel_refit_motion_comparison`** — settles empirically what
the bound can only bracket, by refitting each subject's model three times from
its own `SPM.mat`, changing only the motion block:

- **A** exactly as fitted, **B** the corrected expansion, **C** B without the
  quadratic terms
- everything else held identical, so differences are attributable to the motion
  columns alone
- `r_A_vs_original` validates the rebuild: it must be > 0.999 or nothing else in
  the table means anything

```matlab
OUT = LaBGAScore_firstlevel_refit_motion_comparison(DSGN.modeldir, 'nmotion', 24)
OUT = LaBGAScore_firstlevel_refit_motion_comparison(DSGN.modeldir, 'nmotion', 24, 'variants', {'C'})
```

A and B run by default, because the question the tool exists for is whether an
already-fitted model needs re-running. C is opt-in, and answers the separate
question of which corrected model to re-run *with*; it is ignored when the
existing model has 12 motion columns, where B and C coincide.

`nmotion` describes the **existing** model, not the variants. Datasets are left
clean: unzipped images are removed on exit, the output tree self-ignores, and
the working directory is restored.

## Results per study

### proj_erythritol_4b — 30 subjects, 172 runs, 12 motion columns, HPF 180

| Condition | mean r(task, FD) | t across subjects | exposure |
|---|---|---|---|
| four taste conditions | −0.038 to −0.044 | −5.4 to −7.7 | ~0.09 |
| water | −0.042 | −7.7 | 0.087 |
| `swallow_rinse` | **+0.237** | **+10.6** | 0.246 |

Leakage 0.147–0.174, worst subject 0.26. Consistent but small, and negative:
participants are stiller while the solution is in the mouth, with the movement
parked in `swallow_rinse`, which is modelled as a condition and absorbs it.

### proj_discoverie — MIST, 4 runs, 24 motion columns, HPF 300

| Condition | mean r(task, FD) | t across subjects | exposure |
|---|---|---|---|
| stress | **+0.081** | **+5.86** | **0.619** |
| control | **−0.077** | **−9.41** | **0.505** |

Leakage: stress 0.589, control 0.526, **stress vs control 0.617** — the
differential contrast is the *worst*, because the two conditions' artifacts
have opposite sign and therefore add rather than cancel.

Refit on four subjects, `r_A_vs_original = 1` and `maxdiff = 0` throughout:

| Subject | r(A,B) for stress vs control | median diff | voxels flipping p<.001 |
|---|---|---|---|
| KUL023 | 0.661 | 2.43 SE | 40% |
| CFS0001 | 0.653 | 1.81 SE | 27% |
| CFS0031 | 0.644 | 2.97 SE | 47% |
| CFS0036 | 0.690 | 1.78 SE | 29% |

Implied artifact-to-noise ratio 2.9–4.8, stable across contrasts, so the bound
and the refit agree. The extra collinearity B introduces is confined to the task
regressors — the *nuisance* VIF maxima are comparable in both (144 in A, 129 in
B) — which is the signature of the confound being modelled rather than of a
broken design.

**Two VIF definitions are in play, and their levels are not comparable.**
`scn_spm_design_check` with `'events_only'` scores the task regressors in a
reduced design and gives A 2.09–2.29 against B 7.22–8.76. The refit tool scores
them in the full unfiltered design and gives A 5.2–6.7 against B 13.2–20.6. The
*ratios* agree at about 3× in both. Compare ratios across the two sources, never
absolute levels.

#### Variant C — dropping the quadratics

Three subjects, `stress vs control`:

| Subject | r(B,C) | medabs AB → AC | flip AB → AC | vif_B → vif_C |
|---|---|---|---|---|
| CFS0001 | 0.944 | 1.81 → 1.80 | 0.266 → 0.265 | 17.6 → 13.7 |
| CFS0031 | 0.948 | 2.97 → 2.86 | 0.468 → 0.431 | 19.5 → 12.0 |
| CFS0036 | 0.969 | 1.78 → 1.71 | 0.290 → 0.285 | 13.2 → 11.2 |

C keeps essentially all of the correction — `r_B_vs_C` 0.94–0.97, with
displacement and flip fraction within a few percent of B's — at about
three-quarters of B's task VIF, and with 12 fewer nuisance columns per session,
so ~48 more residual degrees of freedom over four runs. It does **not** return to
A's conditioning: C sits at roughly 2× `vif_A`, because the collinearity is
carried by the linear position terms, which are the whole point of the
correction. The quadratics cost precision without changing the answer.

**C is the model to use for the re-analysis**, and it is what erythritol already
uses (`12hmp`).

Both VIF columns in that table predate the projection-based VIF fix, so the
within-row comparison is internally consistent but the levels are provisional;
`r`, `medabs` and `flip` do not depend on it. The A/B run was repeated afterwards
with the fixed code and gives the values quoted above.

### proj_thc — 13 subjects, 26 runs, two sessions, food images

| Condition | mean r(task, FD) | t across subjects | exposure |
|---|---|---|---|
| thc high / low / neutral | +0.012 / +0.047 / −0.023 | 0.70 / 1.65 / −1.44 | 0.16 / 0.21 / 0.14 |
| pla high / low / neutral | +0.040 / +0.037 / −0.031 | 1.20 / 1.83 / −1.73 | 0.18 / 0.21 / 0.17 |

Leakage 0.214–0.286 across all twelve within-session contrasts, tightly
clustered. Nothing reaches uncorrected significance; smallest p is 0.092.

The Fisher-z between-subject SD is about 0.103 against erythritol's 0.028, so
the coupling is genuinely inconsistent in direction rather than merely
undersampled — the point estimates are the same size as erythritol's, which
reached |t| = 7.7 at n = 30.

Matched thc/placebo pairs share sign in all three cases, so the artifact largely
**cancels** in a drug-versus-placebo difference. Those cross-session contrasts
come back `NaN` from the diagnostic, which works one run at a time; only the
refit can score them.

## Verdicts

| Study | Pattern | Action |
|---|---|---|
| **discoverie** | consistent **and** large; refit moves the primary contrast to r = 0.66 with ~36% of voxels flipping | Refit all subjects, rerun the group analysis. Hold claims about `stress vs control` until then |
| **erythritol_4b** | consistent but small | No re-analysis indicated. Differential contrasts safest |
| **proj_thc** | inconsistent, moderate leakage | Nothing to correct. Expected cost is lost power, not a wrong answer. Cross-session contrasts likely safest of all |

Why discoverie and not the others: MIST has no motion-absorbing nuisance
condition (erythritol's `swallow_rinse` does that job), its long blocks put task
variance at low frequency where the missing subspace also lives, HPF 300 retains
that band, and 24 motion columns make the missing subspace larger.

## What is settled and what is not

**Settled.** The defect is real and its mechanism is understood. The diagnostic
and the refit agree on discoverie. The refit rebuild is faithful
(`r_A_vs_original = 1` in every subject and contrast tested). erythritol and
proj_thc do not need re-analysis. **Variant C is the corrected model to use**:
it reproduces B (`r_B_vs_C` 0.94–0.97) at ~75% of B's task VIF and 12 fewer
nuisance columns per session.

**Not settled.**

- **discoverie's group maps have not been recomputed.** Single-subject con
  images move a great deal; whether the group conclusions change is unknown, and
  should not be guessed either way.
- **The A-versus-B difference conflates bias with variance.** B's stress SE is
  ~1.9× A's, so part of that median 2.4 SE displacement is added noise rather
  than removed bias. C narrows but does not close that gap. The two separate only
  at the group level: bias moves the mean, variance inflates the SE.
- **`vif_C` was measured before the VIF fix** and should be re-read on the next
  run. It changes how much better C is, not whether C wins — that rests on
  `r_B_vs_C`, which does not depend on the VIF computation.
- **Whether the `sub-CFS*` subjects are a discoverie cohort or a separate
  `proj_cfs`** was never confirmed. It changes whether the refit covers one
  study or two.
- **proj_cfs, if separate, has not been assessed at all.**

## Repository state

| Branch | Contents |
|---|---|
| `main` | the two READMEs, every fix except the derivatives one, both audit tools |
| `firstlevel-fixes` | `main` plus one commit: the four `varfun(deriv,…)` → `varfun(zscore,…)` lines, nothing else |

`CANlab_help_examples` `master` carries a pointer to the first-level guides.

The branch is deliberately minimal so it can be reviewed as a single decision.
Merging it changes the noise model for everything fitted afterwards.

Note that merging `firstlevel-fixes` gives model **B**, not C: the fix restores
the parameters alongside their derivatives, and the quadratics are a separate
switch (`LaBGAS_options.movement_reg_quadratic`). Studies following the C verdict
should turn that switch off as well.

## Next steps

1. **Refit discoverie fully** with variant C, and rerun the group analysis.
   Compare group maps, not single-subject con images — that is the only level at
   which bias and added variance separate.
2. **Decide on `firstlevel-fixes`.** Four subjects agreeing removes the reason to
   keep it pending; at minimum it should apply to anything fitted from here on.
3. **Confirm the CFS cohort question**, and run the diagnostic on `proj_cfs` if
   it is a separate study.
4. Only then consider what, if anything, needs communicating externally about
   discoverie.
