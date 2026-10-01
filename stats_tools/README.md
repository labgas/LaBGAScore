# stats_tools — statistics helpers used across the pipeline

Four MATLAB functions that are not specific to any one analysis, plus SAS code
in two flavours: macros for mixed-model effect sizes, and a worked example of
multiple-comparison correction with PROC MULTTEST.

```
functions/              4 MATLAB helpers (below)
sas_macros/             mixed_effectsize.sas, es_identify.sas  -> own README
sas_code/               proc_multtest_fdr.sas, a worked FDR example for the lab
```

**[`sas_macros/README.md`](sas_macros/README.md) is authoritative for the
effect-size macros** — a different language, its own conventions, and its own
caveats.

---

## `functions/`

| function | |
|---|---|
| `LaBGAScore_Storey_FDR.m` | Storey positive-FDR q-values, with a π₀ estimate that is **checked rather than trusted** |
| `LaBGAScore_combat_fit.m` | fits ComBat harmonisation parameters on a **training set only** |
| `LaBGAScore_combat_apply.m` | applies those fitted parameters to held-out data |
| `LaBGAScore_dummy_code.m` | dummy-codes a k-level factor into k−1 indicator columns |

Each carries a full header; what follows is only what you would want to know
before choosing to call one.

### `LaBGAScore_Storey_FDR`

Rewritten because the previous version called `mafdr`, accepted its π₀ unless it
exceeded 0.99, and floored q at p. On real `proj_cfs` data that produced

- 8 ROI SVM p-values → π₀ = 0.895, fine;
- 8 ROI GLM p-values → π₀ = 0.012, i.e. a claim that ~99% of 8 tests are
  non-null, from p-values between .04 and .52.

In the second case every q fell below its own p, so the `q >= p` floor — correct
in itself — handed back the **raw p-values under the name q**. An uncorrected p
reported as an FDR q is the worst failure available here, and it is
**data-dependent**: the same code gave a sane answer on the other set, so one
run cannot reveal it. The old 0.99 guard only ever caught the *conservative*
direction (π₀ → 1, where Storey harmlessly degenerates to BH); the damaging
direction is π₀ → 0.

Four methods: `'sas'` (default), `'lambda'`, `'spline'`, `'bh'`.

**The spline is implemented here, and validated against SAS (2026-10-01).** It
used to delegate to `mafdr`, which is a *different* estimator — and that was the
reason the function disagreed with PROC MULTTEST while labelling its output
"SAS: spline". On two 7-p-value sets:

| p-values | SAS (METHOD=SPLINE) | this function | `mafdr` (old) |
|---|---|---|---|
| `.0063 .0046 .1097 .0017 .0190 .0025 .8031` | 0.17670 | **0.17670** | 0.00319 |
| `.3913 .0124 .2928 .1349 .2515 .0839 .7543` | 0.02547 | **0.02547** | 0.00677 |

n×π₀ reproduces SAS's "estimated number of true nulls" too (1.23693, 0.17829).
Three details decide it, and all three were wrong before: the λ grid is
`(0:19)/20` (SAS's `NLAMBDA=20`), not the `'lambda'` option's default
`0.2:0.1:0.5`; the smoother is a natural cubic **smoothing** spline with 3
effective df; and the estimate is read at the **last λ (0.95), not at λ=1**.
Storey & Tibshirani's paper writes π₀ = s(1), while the qvalue package and SAS
both read off max(λ) — at λ=1 those two sets give 0.14133 and **−0.02118**.
`info.spline_pi0_at_lambda1` reports it so the choice stays visible.

**The `#{p>0.05}` sanity benchmark.** `info.pi0_benchmark` reports
`#{p>0.05}/(n·0.95)`, clamped at 1, and π₀ falling below a third of it earns an
entry in `info.reasons`. This is not a competing estimator — it's a cheap check
on whichever π₀ the chosen method returned.

It works because it *is* Storey's estimator at λ = 0.05, resting on the same
fact: under the null p ~ Uniform(0,1), so only (1−λ) of the true nulls land above
λ. That is what the `/(1−0.05)` divisor corrects for; using `#{p>0.05}/n`
instead undercounts the nulls — a 5% error at λ = 0.05, but a factor of two at
λ = 0.5.

It is a benchmark and **not** a replacement because small p is where the
*alternatives* live: every underpowered true effect is counted as a null, so the
estimate is biased **upward**. That makes it a one-sided reference —

| if π₀ is | then |
|---|---|
| **far below** the benchmark | suspect — and it is the direction that inflates significance |
| above the benchmark | usually benign, merely conservative |

— and a coarse one: at n = 8 it moves in steps of 1/(8·0.95) = 0.13, so it
catches order-of-magnitude disagreement, not fine differences. Note it estimates
a *mixture proportion* and never asserts any individual null is true, which is
why estimating π₀ at all is legitimate frequentist practice.

Calibrated on the six cases where the right answer is known:

| case | n | #p>.05 | benchmark | π₀ before | π₀ after |
|---|---|---|---|---|---|
| SAS example 1 | 7 | 2 | 0.301 | 0.003 | 0.177 |
| SAS example 2 | 7 | 6 | 0.902 | 0.007 | 0.025 |
| proj_cfs m2b | 8 | 7 | 0.921 | 0.012 | 0.000 |
| proj_cfs m2c | 8 | 8 | 1.053 | 0.246 | 0.863 |
| moodbugs roi{2} | 8 | 6 | 0.789 | 0.034 | 0.570 |
| moodbugs roi{4} | 8 | 5 | 0.658 | 0.128 | 0.815 |

All six "before" values sit below a third of their benchmark, so **this check
would have caught the `mafdr` defect on its first run**. After the fix it still
fires on SAS example 2 and proj_cfs m2b — correctly: both have an empty upper
tail and π₀ is not identifiable in either.

**Where π₀ is identified, and when it isn't.** π₀ is read off the *top* of the λ
curve. With no p-value above the last λ, π̂₀(λ) is exactly 0 there and any
estimator extrapolating into that region collapses toward 0 — which is what
produced `mafdr`'s 0.003 above (max p = 0.8031). That condition is now reported
explicitly in `info.reasons` and raised as a warning. Both validated examples
trip it, and SAS's own answer on the second (π₀ = 0.025, i.e. 0.18 of 7
hypotheses null) is not credible either: agreement with SAS is exact, which is
not the same as either number being usable. At small n with a thin upper tail,
BH is the defensible choice.

> **The guard is OFF under the default method.** `'sas'` reproduces SAS PROC
> MULTTEST's PFDR exactly — which is the point, since LaBGAS cross-checks
> analyses against SAS — and that fidelity includes not rejecting a π₀ the
> lambda curve does not support. Under `'sas'` the reliability flag in `info` is
> **advisory: reported, never enforced**. Under the other three methods the
> guard is on and a bad π₀ is rejected. So read the flag, or pick a non-default
> method when SAS agreement is not what you need. Since 2026-10-01 an unreliable
> π₀ that is returned anyway also raises
> `LaBGAScore_Storey_FDR:unreliablePi0`, so it cannot pass into a results table
> silently. The guard is still **not** enforced under `'sas'` — deliberately, so
> the function reproduces PROC MULTTEST including its failures.

### `LaBGAScore_combat_fit` / `_apply`

ComBat split into fit and apply so that harmonisation can happen **inside a
cross-validation fold**: fit on the training fold, apply to held-out subjects.
`combat.m` as distributed fits on everything at once, which leaks held-out
information into the training features.

**Written and documented, but not yet called by anything** — see
[Notes](#notes) for what the decoding script does instead, and why that is
defensible there.

Two properties follow from that purpose and are worth knowing before use:

- **No covariates, deliberately** — `mod = []`, and a non-empty `mod` is
  rejected. This is structural, not a simplification: `combat.m` builds
  `stand_mean` from a design with the batch columns zeroed, so with covariates a
  subject's standardisation depends on that subject's own covariate values, and
  fitted parameters cannot be applied to a held-out subject. It also makes the
  harmonisation **label-blind**, the same invariant the decoding pipelines rely
  on — see [`decoding_toolbox/README.md`](../decoding_toolbox/README.md).
- **Every batch at apply time must have been present at fit time.** A site the
  training fold never saw has no estimated parameters, so **leave-one-site-out
  designs are not supported**. The R reference implementation raises the same
  error.

### `LaBGAScore_dummy_code`

A k-level factor needs k−1 columns. Coding three sites as a single −1/0/1 column
— as `proj_discoverie`'s `model_3a` does for `center` — treats them as *ordered*
and spends one degree of freedom, so it removes only part of the between-site
variance and imposes an ordering that does not exist. Any nuisance adjustment
for an unordered factor should go through this function.

Accepts numeric, logical, char, cellstr, string or categorical; returns double
indicator columns plus the level names and the input as a categorical.

## `sas_code/proc_multtest_fdr.sas` — and the one check to do every time

A worked example running five real p-value sets through

```sas
proc multtest inpvalues=<data> fdr pfdr afdr plots=all;
```

which puts three corrections side by side — **FDR** (Benjamini–Hochberg, assumes
π₀ = 1), **PFDR** (Storey, scales BH by an *estimated* number of true nulls m₀),
**AFDR** (adaptive, estimates m₀ by LOWESTSLOPE (Benjamini & Hochberg 2000), built for small m) — plus the
diagnostic plots, including the fitted m₀.

**Everything turns on m₀, because q = (m₀/m) · q_BH.** Halve m₀ and every q
halves. So check it, using the benchmark described above:

> **How many p-values exceed 0.05?** That count, over 0.95, is itself a rough
> estimate of m₀, biased upward. An m₀ far *below* it is the direction that
> invents significance; at or above it is merely conservative.

Two estimates are wrong on their face regardless: **m₀ = 0** asserts no
hypothesis is null, and **m₀ < 1 with m < 10** claims less than one null among a
handful of tests.

### When the check fails

Syntax below is verified against the
[SAS/STAT 14.1 MULTTEST documentation](https://support.sas.com/documentation/onlinedoc/stat/141/multtest.pdf).

1. **Try another m₀ estimator.** The option is `NTRUENULL=`, with **`M0=` as a
   documented alias**; `PTRUENULL=` takes a *proportion* instead of a count, and
   either accepts a plain number in place of a keyword. Keywords: `SPLINE`,
   `BOOTSTRAP`, `LOWESTSLOPE`, `DECREASESLOPE`, `LEASTSQUARES`, `KSTEST`,
   `MEANDIFF`. The defaults differ **by adjustment**, which is worth knowing
   before concluding two adjustments disagree about the data:

   | adjustment | default m₀ estimator |
   |---|---|
   | `PFDR` | `SPLINE`, falling back to `BOOTSTRAP` if the estimate is nonpositive or the spline's slope at the last λ exceeds 0.1 × the range of fitted values |
   | `ADAPTIVEFDR` (alias `AFDR`) | `LOWESTSLOPE` |
   | `ADAPTIVEHOLM`, `ADAPTIVEHOCHBERG` | `DECREASESLOPE` |

   `LOWESTSLOPE` and `DECREASESLOPE` are the ones to reach for at small m: they
   read m₀ off the slope of the ordered p-values and need no well-populated upper
   tail — exactly what `SPLINE` and `BOOTSTRAP` lack when only a handful of
   p-values exceed 0.05.
2. **Or report `ADAPTIVEFDR`**, but only if *its* m₀ passes the same check. It
   defaults to `LOWESTSLOPE` for this reason, and in all five examples it passes
   whenever PFDR fails.
3. **Or report plain `FDR`.** Always defensible, and the honest answer when m₀ is
   not estimable — which it is not when only a handful of p-values exceed 0.05.

### Reading the diagnostic plots

`PLOTS=ALL` includes the two that matter:

| plot | what it is |
|---|---|
| **LambdaPlot** | *"MSE or NTRUENULL by lambda"* — produced for `PFDR`, or `NTRUENULL=SPLINE`/`BOOTSTRAP`. **This is the λ curve the estimate is read off**, so look at its right-hand end, where m₀ is taken: flat and settled means identifiable; erratic, rising, or implying a proportion above 1 means it is not, and no choice of λ rescues it |
| **RawUniformPlot** | raw p-values by rank plus their histogram. Near-uniform under a mostly-null family. Almost nothing above 0.05 (see K1ROI) says a small m₀ is plausible; flat across [0,1] says m₀ should be near m |

Do **not** set m₀ by hand without external grounds for the number: m₀ = m is
exactly BH, and anything lower buys significance by assumption rather than from
the data.

### What the five examples show

Computed with `LaBGAScore_Storey_FDR`, which reproduces PROC MULTTEST's PFDR
spline exactly (verified on the first two to five decimals):

| dataset | m | #p>.05 | benchmark m₀ | PFDR m₀ | AFDR m₀ | verdict |
|---|---|---|---|---|---|---|
| CytokinesT1 | 7 | 2 | 2.1 | 1.24 | 6.00 | both pass |
| CytokinesT2 | 7 | 6 | 6.3 | **0.18** | 5.00 | PFDR fails → report AFDR |
| VTROI | 14 | 12 | 12.6 | 7.50 | 13.00 | both pass |
| K1ROI | 14 | 1 | 1.1 | **0.00** | 14.00 | PFDR m₀ = 0, impossible |
| SCFAs | 4 | 4 | 4.2 | **1.06** | 4.00 | PFDR fails → report FDR |

**CytokinesT2** is the clearest failure: six of seven p-values exceed 0.05, yet
PFDR puts the number of true nulls at 0.18. AFDR's 5 is credible.

**SCFAs** shows the floor of the method — m = 4, all above 0.17, and PFDR claims
1.06 nulls. Four tests cannot support estimating m₀; report FDR.

**K1ROI** is the instructive edge case. Thirteen of fourteen p-values are *below*
0.05, so the benchmark is low (1.1) and a small m₀ is genuinely plausible — but
PFDR returns exactly 0, which cannot be right, while AFDR's 14 is
over-conservative given the data. When the two bracket the answer that widely,
say which you used and why.

## Dependencies

No CanlabCore or SPM dependency. `LaBGAScore_Storey_FDR` uses MATLAB's
Statistics and Machine Learning Toolbox for the spline fit. The ComBat functions
mirror the arithmetic of `combat.m` (Johnson, Li & Rabinovic 2007; Fortin's
MATLAB port) but do not call it.

## Notes

- These are **library functions** — safe to fix in place, and a behaviour change
  reaches every caller in every repo. Current callers:

  | function | called from |
  |---|---|
  | `LaBGAScore_Storey_FDR` | `decoding_toolbox/LaBGAScore_decoding_SVM_between_subjects.m`; `prep_3a_run_second_level_regression_and_save.m` and `h_signature_responses_group_diff.m` in the fork |
  | `LaBGAScore_dummy_code` | the same decoding script, and `prep_3a` in the fork |
  | `LaBGAScore_combat_fit` / `_apply` | **no caller yet** — see below |

- **The ComBat fit/apply pair is implemented but not yet wired in.** The
  decoding script names both in its comments and then calls plain `combat(...)`
  directly, harmonising once over all subjects rather than per fold. That is a
  defensible choice there and its header argues it: harmonisation is label-blind
  (`mod = []`) and applied identically to the real and permuted runs, so it
  cannot leak *labels* even when fitted on everything. The fit/apply pair exists
  for the stricter regime — where you want no held-out subject to influence the
  training-fold features at all — and is ready for a caller that wants it.
- `LaBGAScore_dep_map.m` walks `.m` files only, so `sas_macros/` is covered by
  no dependency tooling and appears in no generated documentation.
- The SAS macros' own README records that they have **not yet been executed
  against a real model**: their arithmetic was verified independently and their
  syntax checked by reading, no more. Remove that caveat there once someone has
  run them.
