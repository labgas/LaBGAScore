# stats_tools — statistics helpers used across the pipeline

Four MATLAB functions that are not specific to any one analysis, plus the SAS
macros for mixed-model effect sizes.

**[`sas_macros/README.md`](sas_macros/README.md) is authoritative for the SAS
side** (`mixed_effectsize.sas`, `es_identify.sas`) — a different language, its
own conventions, and its own caveats. This file covers `functions/` only.

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

> **The guard is OFF under the default method.** `'sas'` reproduces SAS PROC
> MULTTEST's PFDR exactly — which is the point, since LaBGAS cross-checks
> analyses against SAS — and that fidelity includes not rejecting a π₀ the
> lambda curve does not support. Under `'sas'` the reliability flag in `info` is
> **advisory: reported, never enforced**. Under the other three methods the
> guard is on and a bad π₀ is rejected. So read the flag, or pick a non-default
> method when SAS agreement is not what you need.

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
