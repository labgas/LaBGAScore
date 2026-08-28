# sas_macros

SAS macros for statistics that the MATLAB side of LaBGAScore does not cover.
Currently: effect sizes for fixed effects in models fitted with `PROC MIXED`.

| File | Macros | Purpose |
|---|---|---|
| `mixed_effectsize.sas` | `%mixed_effectsize` | Effect sizes for a fitted model |
| `es_identify.sas` | `%es_identify`, `%es_identify_ds` | Audit effect sizes in an existing results table |

Both files are self-contained. `%include` the one you need; there are no
cross-dependencies.

> **Status: not yet run in SAS.** The arithmetic has been verified
> independently and the macro syntax statically checked, but neither file has
> been executed against a real model. Run them against a model with known
> output and confirm the results before using them in a manuscript, then
> delete this notice.

---

## Why these exist

An earlier in-house script (*Calculation of effect size following marginal
linear mixed models*, B. Dalile, 3 December 2021) computed three statistics —
`eta_2`, `omega_2` and `partial_eta_2` — and it was easy to report one of them
under the name of another. It also contained a bug. The error sum of squares
was formed as

```sas
ss_error = mse * (_freq_ - numdf);
```

which uses the **number of observations** minus the effect's **numerator**
degrees of freedom, where the model's **denominator** degrees of freedom
belong. Substituting the definitions, the `mse` cancels out entirely and the
expression reduces to

```
numdf*F / (numdf*F + N_obs - numdf)
```

In a repeated-measures design with ~660 observations and a denominator df of
~161, that inflates the denominator roughly fourfold and shrinks partial
eta-squared by the same factor.

Partial eta-squared needs **nothing but the Type 3 table**, because `PROC
MIXED` already reports the Kenward-Roger or Satterthwaite denominator df
there:

```
partial_eta2 = (NumDF * FValue) / (NumDF * FValue + DenDF)
```

Neither macro here reads `_FREQ_`.

---

## `mixed_effectsize.sas`

### Usage

```sas
%include "path/to/sas_macros/mixed_effectsize.sas";

proc mixed data = mydata;
    class id group time;
    model y = group|time / outpm = outpm_y ddfm = kr;
    repeated time / subject = id type = un;
    ods output tests3 = tests3_y;
run;

%mixed_effectsize(tests3 = tests3_y,
                  model  = MARGINAL,
                  resids = outpm_y,
                  dv     = y,
                  out    = es_y,
                  label  = %str(Self-reported stress, MAST));
```

Partial eta-squared and its confidence interval need only `tests3`. Supply
`resids=` and `dv=` as well if eta-squared and omega-squared are also wanted.

### Output

Per effect: `partial_eta2` with a noncentral-F confidence interval,
`cohen_f2`, and — only when a residuals dataset is supplied — eta-squared and
omega-squared, labelled as approximations. Every quantity has its own variable
name and a full printed label, which is the point.

### Marginal versus mixed models

`model = MARGINAL | MIXED`.

**Partial eta-squared is identical either way.** The formula uses only
`NumDF`, `DenDF` and `FValue`, and the denominator df already reflects
whatever covariance structure was fitted. Nothing to choose.

What changes is eta-squared, because "the residual" means two different
things:

| `residtype` | Residual | `SS_error` contains | eta-squared is | Output name |
|---|---|---|---|---|
| `MARGINAL` (default) | y − Xb, from `OUTPM=` | between-subject variance | proportion of **total** variance | `eta2_approx` |
| `CONDITIONAL` | y − Xb − Zu, from `OUTP=` | within-subject variance only | proportion of **within-subject** variance | `eta2_cond_approx` |

Conditional eta-squared is systematically larger and is **not comparable with
eta-squared from other studies**, because the random effects have already
absorbed the between-subject variance. The macro therefore defaults to
marginal residuals for both model types, names the two outputs differently so
they cannot be confused in a results table, and refuses
`residtype = CONDITIONAL` when `model = MARGINAL` — a model with no `RANDOM`
statement has no BLUPs, so `OUTP=` and `OUTPM=` are identical.

The macro cannot tell `OUTPM=` from `OUTP=` by inspection, so it prints a NOTE
stating which it assumed rather than guessing. Read the log.

### Which statistic to report

Report `partial_eta2`. It is the convention in this literature, and a reader
can recompute it from the F and df printed in the same table — which is what
makes a results table checkable.

Eta-squared cannot be recovered from F and df alone, so a reader cannot verify
it. Report it only alongside partial eta-squared, never instead of it, and
name it correctly.

For a two-group contrast specifically, a standardised mean difference with a
confidence interval (`ESTIMATE` or `LSMESTIMATE` with a sensible divisor) is
more directly interpretable than any variance-explained measure.

### Caveats to state in a Methods section

1. Eta-squared and omega-squared are **not uniquely defined** for a mixed or
   marginal model — there is no single residual variance and no single error
   df. They are formed here as `SS_error = USS(residuals)`,
   `SS_total = CSS(dv)`, `MSE = SS_error / DenDF`,
   `SS_effect = NumDF * FValue * MSE`. That is the natural analogue of the
   fixed-effects definitions but remains an approximation, hence the
   `*_APPROX` names.
2. Under Kenward-Roger, `DenDF` differs between effects, so the implied MSE
   differs slightly too. Partial eta-squared is unaffected (the MSE cancels);
   eta-squared and omega-squared are mildly affected. The macro reports the
   `DenDF` used for each row so the dependence is visible.
3. With `residtype = CONDITIONAL`, report the variance components alongside
   the effect sizes. Without them a reader cannot tell how much variance the
   random effects absorbed.
4. The confidence interval inverts the noncentral F distribution
   (Steiger 2004) using N = `NumDF + DenDF + 1`. With Kenward-Roger df that N
   is an approximation, so the interval is approximate.
5. **Use `DDFM = KR`** (or `KR2`). With the default containment method the
   denominator df can be badly wrong for these designs, and every quantity
   here depends on it. The macro warns if every `DenDF` is a whole number,
   which usually means KR or Satterthwaite was not requested.

### Correspondence with `PROC GLM`

`PROC GLM`'s `EFFECTSIZE` option defines the **partial correlation ratio** as
`SS_effect / (SS_effect + SS_error)`. The formula used here is algebraically
identical. In GLM, `F = (SS_effect/df_effect) / MSE` with
`MSE = SS_error/df_error`, so `SS_effect = df_effect·F·MSE` and
`SS_error = df_error·MSE`. Substituting, the MSE cancels:

```
SS_effect/(SS_effect + SS_error) = (df_effect·F·MSE) / (df_effect·F·MSE + df_error·MSE)
                                 = (NumDF·F) / (NumDF·F + DenDF)
```

Verified numerically on an unbalanced 2×3 factorial with Type III sums of
squares:

| Effect | NumDF | F | SAS definition | Formula used here | Difference |
|---|---|---|---|---|---|
| group | 1 | 0.9276 | 0.01935469684 | 0.01935469684 | 0 |
| time | 2 | 10.8113 | 0.31509459617 | 0.31509459617 | 5.6e-17 |
| group×time | 2 | 7.8300 | 0.24992026401 | 0.24992026401 | 2.8e-17 |

Semipartial (plain) eta-squared, `SS_effect/SS_total`, likewise agrees exactly
— confirming that reconstructing `SS_effect` as `NumDF·FValue·MSE` is exact in
the fixed-effects case.

**Two things "identical formula" does not mean:**

1. In GLM, `DenDF` is the single residual df, the same for every effect. Under
   Kenward-Roger in `PROC MIXED` it is effect-specific and usually fractional.
   That is correct — it is the right df for that effect's test — but the values
   are not slices of one sum-of-squares decomposition the way GLM's are.
2. In a mixed or marginal model there is **no exact sum-of-squares
   decomposition at all**. With correlated errors the F test is a Wald-type
   test on generalised least squares estimates, not a ratio of orthogonal sums
   of squares. The identity therefore runs one way only: the formula can be
   computed, but the result cannot be recovered as
   `SS_effect/(SS_effect + SS_error)` from any real decomposition. Read it as
   *the effect size implied by the F test* — which is precisely what Edwards et
   al. (2008) formalise as R²_β, and is worth one sentence in a Methods
   section.

### Open question: which omega-squared convention

The omega-squared computed here uses

```
(SS_effect − NumDF·MSE) / (SS_total + MSE)
```

inherited from the 2021 script. That is the classic Hays/Keppel convention, and
is what Olejnik & Algina (2003) use.

The SAS documentation for `EFFECTSIZE` renders semipartial omega-squared with
`SS_total` alone in the denominator in one place, and describes it as "the total
sum of squares adjusted by MSE terms" in another. **The formulas on those pages
are published as images rather than text, so this has not been verified against
SAS's actual output.**

The difference is small — about 0.8% in a test case — and partial eta-squared is
unaffected either way. To settle it, run `PROC GLM` with `EFFECTSIZE` on any
small dataset and compare the printed semipartial omega-squared against both

```
(SS_effect − df_effect·MSE) / SS_total
(SS_effect − df_effect·MSE) / (SS_total + MSE)
```

then update this README and the macro header with what you find.

### References

- Edwards LJ, Muller KE, Wolfinger RD, Qaqish BF, Schabenberger O. An R² statistic
  for fixed effects in the linear mixed model. *Statistics in Medicine*
  2008;27(29):6137–6157. [doi:10.1002/sim.3429](https://doi.org/10.1002/sim.3429)
  (PMID 18816511)
- Olejnik S, Algina J. Generalized eta and omega squared statistics: measures of
  effect size for some common research designs. *Psychological Methods*
  2003;8(4):434–447.
  [doi:10.1037/1082-989X.8.4.434](https://doi.org/10.1037/1082-989X.8.4.434)
- Steiger JH. Beyond the F test: effect size confidence intervals and tests of
  close fit in the analysis of variance and contrast analysis. *Psychological
  Methods* 2004;9(2):164–182.
  [doi:10.1037/1082-989X.9.2.164](https://doi.org/10.1037/1082-989X.9.2.164)
  (PMID 15137887)
- SAS Institute. Effect Size Measures for F Tests in GLM. *SAS/STAT User's
  Guide*.
  [support.sas.com](https://support.sas.com/documentation/cdl/en/statug/63962/HTML/default/statug_glm_sect032.htm)

---

## `es_identify.sas`

Forensic helper for auditing effect sizes in tables produced by other code,
including the older script described above.

### The decisive test

For an effect with numerator df *q*, F statistic *F*, and reported effect size
*r*: if the reported values are eta-squared, then

```
r / (q * F)
```

must be **constant across every effect within one model**, whatever the sample
size happens to be. This requires no assumptions at all — not even N.

Agreement to within a few percent (allowing for rounding in a published table)
means the reported values are eta-squared. Disagreement by a factor of two or
more means the value is something else and must be checked against the
original output. An implied residual/total variance ratio above 1 is
impossible and is conclusive.

### Usage

One effect at a time, printed to the log:

```sas
%es_identify(numdf = 1, f = 40.45, dendf = 161, nobs = 664, reported = 0.028);
```

A whole model at once, which is how the test is meant to be used:

```sas
data mast_stress;
    length effect $32;
    infile datalines dsd;
    input effect $ numdf dendf f reported;
    datalines;
group,1,161,40.45,0.028
time,3,161,66.22,0.130
group*time,3,161,2.44,0.005
;
run;

%es_identify_ds(data = mast_stress, nobs = 664);
```

`nobs=` is optional; without it the macro still reports the correct partial
eta-squared and the sample-size-free ratio, but cannot evaluate the old
script's formula or the implied variance ratio.

---

## Conventions

- Two-space indentation inside macros, aligned `=` in parameter lists.
- Every macro validates its inputs and aborts with an `ERROR:` in the log
  rather than producing a silently wrong dataset.
- Every macro states its assumptions in the log with `NOTE:` when it cannot
  verify them from the data — in particular, which residual type it assumed.
- Approximations are named `*_APPROX`. Anything that cannot be recomputed by a
  reader from the printed test statistics is flagged as such in the header.
