# Motion regressor audit — state of play

Working note, 2026/09/03. Written as a handover: enough to pick this up cold, on
another machine, without the conversation it came from.

The short version: the first-level motion regressors have never been what the
script headers say they are. Whether that matters turns out to depend entirely
on the study, and of the three checked so far, one needs a re-analysis, one does
not, and one is unaffected.

For the study that does — discoverie — the single-subject con images move a
great deal, and on the primary contrast that movement is **removed bias rather
than added noise**. The model to re-run it with is **variant C**: the corrected
motion expansion without the quadratic terms. What has not been done is the
group analysis, which is the only thing that can say whether the published
conclusions change.

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
the bound can only bracket, by refitting each subject's model from its own
`SPM.mat`, changing only the motion block:

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

Three subjects, all VIFs from the projection-based computation:

| Subject | contrast | r(A,B) | r(A,C) | r(B,C) | medabs AB → AC | vif_A | vif_B | vif_C |
|---|---|---|---|---|---|---|---|---|
| CFS0001 | stress vs control | 0.653 | 0.655 | 0.944 | 1.81 → 1.80 | 5.32 | 17.66 | 13.35 |
| CFS0031 | stress vs control | 0.644 | 0.696 | 0.948 | 2.97 → 2.86 | 6.72 | 20.62 | 11.92 |
| CFS0036 | stress vs control | 0.690 | 0.713 | 0.969 | 1.78 → 1.71 | 5.20 | 13.20 | 11.19 |

C keeps essentially all of the correction — `r_B_vs_C` 0.86–0.98 across all nine
rows — at about three-quarters of B's task VIF, and with 12 fewer nuisance
columns per session, so ~48 more residual degrees of freedom over four runs. It
does not return to A's conditioning on `stress` (~2× `vif_A`), because the
collinearity is carried by the linear position terms, which are the whole point
of the correction; on `control` it comes close (1.1–1.6× `vif_A`).

**C is the model to use for the re-analysis**, and it is what erythritol already
uses (`12hmp`).

#### The primary contrast's displacement is bias, not added variance

The A→C displacement is smaller than A→B throughout, which is what dropping 12
collinear columns should do. But *how much* smaller answers the question this
audit has so far been unable to settle — how much of the A-to-B displacement is
removed bias and how much is added noise:

| Contrast | medabs A→B | medabs A→C | shrinkage |
|---|---|---|---|
| `stress` | 1.62 / 3.46 / 2.00 | 1.59 / 3.12 / 1.76 | 2–12% |
| `control` | 1.04 / 1.85 / 1.33 | 0.91 / 1.62 / 1.05 | 12–21% |
| **`stress vs control`** | 1.81 / 2.97 / 1.78 | 1.80 / 2.86 / 1.71 | **0.5–4%** |

C carries a third less collinearity than B, so any part of the A→B displacement
that was inflated noise should shrink with it. On the simple effects it does, by
10–20%. On the **primary differential contrast it does not** — 0.5–4%, and
`r_A_vs_C` (0.655/0.696/0.713) is no higher than `r_A_vs_B` (0.653/0.644/0.690).

The quadratics' added noise is largely common to the two conditions and cancels
in the difference, while the artifact has opposite sign in stress and control and
therefore adds. So for `stress vs control` the con images move because the
artifact is being removed, not because the corrected model is noisier. That is
the contrast the study's conclusions rest on, and it is the one where the
displacement is real.

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
| **discoverie** | consistent **and** large; refit moves the primary contrast to r = 0.66 with ~36% of voxels flipping, and that movement is bias rather than added noise | Refit all subjects with **variant C**, rerun the group analysis. Hold claims about `stress vs control` until the group maps exist |
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
it reproduces B (`r_B_vs_C` 0.86–0.98) at ~75% of B's task VIF and 12 fewer
nuisance columns per session. And for the **primary contrast the displacement is
bias rather than added variance** — it survives the drop in collinearity from B
to C essentially unchanged (0.5–4%), where the simple effects shrink 10–20%.

**Not settled.**

- **discoverie's group maps have not been recomputed.** Single-subject con
  images move a great deal; whether the group conclusions change is unknown, and
  should not be guessed either way.
- **The bias/variance split is resolved for `stress vs control` only.** On the
  simple effects a tenth to a fifth of the A→B displacement is inflated noise,
  and the residue is bounded from one direction only: B's stress SE is ~1.9× A's,
  C's is lower but still above A's. Only the group level separates them fully —
  bias moves the mean, variance inflates the SE.
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
