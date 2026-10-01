/*****************************************************************************
  proc_multtest_fdr.sas

  Multiple-comparison correction with PROC MULTTEST, as a worked example for the
  lab. Each block below feeds one set of raw p-values to

      proc multtest inpvalues=<data> fdr pfdr afdr m0=decreaseslope plots=all;

  which requests three corrections side by side, plus the diagnostic graphics:

      FDR    Benjamini-Hochberg. Assumes every hypothesis could be null
             (pi0 = 1). Always valid, never the most powerful.
      PFDR   Storey's positive FDR. Scales BH by an ESTIMATED number of true
             nulls (m0). More powerful when that estimate is good, and
             anti-conservative when it is not.
      AFDR   Adaptive FDR. Also estimates m0, but applies it as a step-up
             procedure on m0 rather than as Storey's scalar rescaling.
      M0=DECREASESLOPE   the estimator used for m0. THIS IS THE LAB DEFAULT and
             is set deliberately here - see the next section. It overrides the
             per-adjustment defaults, so PFDR and AFDR both use it and their
             m0 agrees; their q-values then usually coincide too.
      PLOTS=ALL  the diagnostic plots, including the one that matters here: the
             fitted estimate of the number of true nulls.

  ---------------------------------------------------------------------------
  WHY M0=DECREASESLOPE RATHER THAN THE SAS DEFAULT
  ---------------------------------------------------------------------------

  Left to itself, PFDR estimates m0 with SPLINE, falling back to BOOTSTRAP. Both
  read m0 off the UPPER tail of the p-value distribution, and at the sample sizes
  this lab works at that tail is nearly empty - so both extrapolate into a region
  with no data, and extrapolate DOWNWARD, which is the direction that invents
  significance.

  Measured, not assumed. 400 simulated datasets per cell, LaBGAScore_Storey_FDR
  (which reproduces PROC MULTTEST's spline exactly):

      m = 8, all nulls, nominal FDR 0.05
          SPLINE default   realised FDR 0.178   degenerate m0 in 25.5% of runs
          DECREASESLOPE    realised FDR 0.065   degenerate m0 in  0.0% of runs

  The spline default does not control FDR at small m. DECREASESLOPE does, and
  still beats plain BH on power (0.711 against 0.674) where there is signal.
  On all five datasets below it also passes the sanity check that follows, which
  the spline fails on three of them.

  So: specify M0=DECREASESLOPE. The same change has been made to the lab's MATLAB
  implementation (LaBGAScore_Storey_FDR, default 'decreaseslope' since
  2026-10-01), so SAS and MATLAB agree out of the box.

  ONE COST, worth knowing: the LambdaPlot described below is the diagnostic for
  the SPLINE/BOOTSTRAP lambda curve, and per the SAS documentation it is produced
  for PFDR's own default or for NTRUENULL=SPLINE|BOOTSTRAP - there is no lambda
  curve to draw for DECREASESLOPE. If you want to SEE why the spline fails on
  your data, run the same block once without M0= and look at the plot; report the
  M0=DECREASESLOPE numbers. (Stated from the documentation, not measured here -
  check your own output.)

  ---------------------------------------------------------------------------
  M0=DECREASESLOPE DOES NOT RETIRE THE CHECK BELOW - IT MAKES IT QUIETER
  ---------------------------------------------------------------------------

  Setting the better estimator is not the same as not having to look. DECREASESLOPE
  fixes the DEGENERATE failures - m0 at or near 0, asserting that essentially
  nothing is null - thoroughly. It does NOT guarantee m0 is well estimated, and at
  these sample sizes it often is not.

  On the five datasets below, DECREASESLOPE's m0 still falls BELOW the #{p>0.05}
  benchmark on THREE of the five (CytokinesT1 2.00 vs 2.1, VTROI 11.00 vs 12.6,
  CytokinesT2 3.00 vs 6.3) - even though none of them is absurd any more. Measured
  more broadly with the MATLAB implementation: m0 below benchmark in 32%% of the 73
  p-value families of the lab's eight final second-level models, and in 41.7%% of
  random panels - against 85.3%% for the spline. The catastrophic cases are gone;
  the mildly anti-conservative ones are not.

  So run the check below on every PFDR or AFDR output, whatever M0= says.

  ---------------------------------------------------------------------------
  AND COMPARE ESTIMATORS - IT COSTS ONE EXTRA PROC STEP
  ---------------------------------------------------------------------------

  One estimator is one opinion. Add a second run with a different NTRUENULL= and
  see whether they agree:

      proc multtest inpvalues=<data> pfdr m0=lowestslope;   run;

  m0 on the five datasets below, by method:

      dataset        m   benchmark   DECREASESLOPE   LOWESTSLOPE   SPLINE
      CytokinesT1    7      2.1          2.00           6.00        1.24
      CytokinesT2    7      6.3          3.00           5.00        0.18
      VTROI         14     12.6         11.00          13.00        7.50
      K1ROI         14      1.1         14.00          14.00        0.00
      SCFAs          4      4.2          4.00           4.00        1.06

  LOWESTSLOPE is the more conservative of the two slope estimators on every one
  of these, and lands ABOVE the benchmark where DECREASESLOPE lands below it.
  Where the two disagree materially - CytokinesT1, 2 against 6 - report which you
  used and why, or fall back to plain FDR. Where they agree, the q-values rest on
  more than a single number. (Remember that specifying M0= applies it to PFDR and
  AFDR BOTH, so a second opinion needs a second PROC step, not a second
  adjustment in the same one.)

  ---------------------------------------------------------------------------
  THE ONE CHECK YOU MUST DO: IS THE ESTIMATED NUMBER OF TRUE NULLS SANE?
  ---------------------------------------------------------------------------

  PFDR and AFDR both report an "Estimated number of true null hypotheses" (m0).
  That number drives every q-value: q = (m0/m) * q_BH. Halve m0 and you halve
  every q. It is therefore the single quantity worth checking, and the check is
  arithmetic you can do in your head:

      BENCHMARK:  how many p-values exceed 0.05?

  Under the null, p is Uniform(0,1), so about 95% of true nulls land above 0.05.
  #{p > 0.05}/0.95 is thus itself a rough estimate of m0 - and it is biased
  UPWARD, because every underpowered real effect (p = .08, p = .12) is counted as
  a null. That makes it a soft UPPER reference:

      m0 far BELOW #{p>0.05}/0.95   ->  DO NOT TRUST IT. This is the direction
                                        that invents significance.
      m0 at or above the benchmark  ->  fine, merely conservative.

  Two estimates are wrong on their face whatever the benchmark says:
      m0 = 0            asserts that NO hypothesis is null. Cannot be right.
      m0 < 1 with m < 10  claims less than one null among a handful of tests.

  ---------------------------------------------------------------------------
  WHAT TO DO WHEN THE CHECK FAILS
  ---------------------------------------------------------------------------

  Syntax below verified against the SAS/STAT 14.1 MULTTEST documentation
  (support.sas.com/documentation/onlinedoc/stat/141/multtest.pdf).

  1. Try a different m0 estimator - the blocks below already do this, via
     M0=DECREASESLOPE. The option is NTRUENULL=, with M0= as a documented alias; PTRUENULL= takes a PROPORTION instead of a count, and
     either accepts a positive integer (resp. a proportion) in place of a
     keyword. The keywords are:

         SPLINE          Storey & Tibshirani (2003) cubic spline. pi0_hat(lambda)
                         = #{p>lambda}/(m(1-lambda)) on lambda in {0,1/n,...,
                         (n-1)/n}, fitted with a natural cubic spline of 3 df
                         (DF= to change), read off at the LAST lambda.
         BOOTSTRAP       Storey & Tibshirani (2003) bootstrap. Picks the lambda
                         minimising MSE. NBOOT=10000 and NLAMBDA=20 by default.
         LOWESTSLOPE     Benjamini & Hochberg (2000).
         DECREASESLOPE   Schweder & Spjotvoll (1982) as modified by Hochberg &
                         Benjamini (1990).
         LEASTSQUARES    least-squares search for the cutpoint.
         KSTEST          Kolmogorov-Smirnov uniformity test, Turkheimer et al.
                         (2001).
         MEANDIFF        mean of differences, Hsueh, Chen & Kodell (2003).

     The defaults differ by adjustment, which is worth knowing before you
     conclude that two adjustments disagree about the data - and worth knowing
     because specifying M0= (as these blocks do) replaces ALL of them at once,
     which is why PFDR and AFDR below report the same m0:

         PFDR            SPLINE first; if the estimate is nonpositive, or if the
                         slope of the spline at the last lambda exceeds 0.1 times
                         the range of the fitted spline values, BOOTSTRAP instead
         ADAPTIVEFDR     LOWESTSLOPE
         ADAPTIVEHOLM    DECREASESLOPE
         ADAPTIVEHOCHBERG  DECREASESLOPE

     LOWESTSLOPE and DECREASESLOPE are the ones to reach for at small m: they
     read m0 off the slope of the ordered p-values and need no well-populated
     upper tail, which is exactly what SPLINE and BOOTSTRAP lack when only a
     handful of p-values exceed 0.05. DECREASESLOPE is the lab default for that
     reason; LOWESTSLOPE is the natural second thing to try, and is the more
     conservative of the two on every dataset below.

  2. Or report ADAPTIVEFDR (alias AFDR) instead - but only if ITS m0 passes the
     same check. Note that once you specify M0=, both adjustments use the
     estimator you named, so this is no longer a way to get a SECOND opinion on
     m0; drop the M0= option, or name a different estimator, to get one.

  3. Or report plain FDR (Benjamini-Hochberg). Always defensible, and the honest
     answer when m0 simply is not estimable - which it is not when only a handful
     of p-values sit above 0.05.

  Do NOT respond to a failed check by specifying m0 by hand unless you have
  external grounds for the number. Specifying m0 = m (or PTRUENULL=1) is exactly
  BH; anything lower buys significance by assumption rather than from the data.

  ---------------------------------------------------------------------------
  READING THE DIAGNOSTIC PLOTS
  ---------------------------------------------------------------------------

  PLOTS=ALL includes the two that matter here:

      LambdaPlot       "MSE or NTRUENULL by lambda". Produced for PFDR, or for
                       NTRUENULL=SPLINE / BOOTSTRAP. This IS the lambda curve the
                       estimate is read off, so look at its right-hand end, where
                       m0 is taken. A curve that is flat and settled there is
                       identifiable; one that is erratic, or rises, or implies a
                       proportion above 1, is not - and then no choice of lambda
                       rescues it.
      RawUniformPlot   raw p-values by rank, plus their histogram. Under a mostly
                       null family this is close to uniform. A histogram with
                       almost nothing above 0.05 (see K1ROI below) tells you a
                       small m0 is plausible; one that is flat across [0,1] tells
                       you m0 should be near m.

  ---------------------------------------------------------------------------
  WHAT TO EXPECT FROM THE FIVE EXAMPLES BELOW
  ---------------------------------------------------------------------------

  Computed with LaBGAScore_Storey_FDR, which reproduces PROC MULTTEST exactly -
  the spline on CytokinesT1/T2 to five decimals, and DECREASESLOPE on both of
  those against SAS output (m0 = 2 and m0 = 3, q-values agreeing to 5e-5):

                                       m0 if you ask for
    dataset        m   #p>.05  bench   DECREASESLOPE  SPLINE   LOWESTSLOPE
    CytokinesT1    7      2      2.1        2.00       1.24       6.00
    CytokinesT2    7      6      6.3        3.00       0.18       5.00
    VTROI         14     12     12.6       11.00       7.50      13.00
    K1ROI         14      1      1.1       14.00       0.00      14.00
    SCFAs          4      4      4.2        4.00       1.06       4.00

  Read the DECREASESLOPE column against the benchmark: it is at or near it on
  four of five, where the spline sits far below on three. That is the whole case
  for the changed default, visible on five real datasets.

  CytokinesT1 - both estimators pass, but they still disagree about the data.
  DECREASESLOPE's m0 = 2 lands exactly on the benchmark. Note that the smallest
  p-value, 0.0017, comes back as q = 0.0025: even a passing m0 moves things.

  CytokinesT2 is the clearest spline failure: six of seven p-values exceed 0.05,
  yet SPLINE puts the number of true nulls at 0.18. DECREASESLOPE says 3. That is
  still BELOW the benchmark of 6.3, so this dataset does not fully pass even
  under the new default - it is merely no longer absurd. Treat the q-values here
  as optimistic and consider reporting plain FDR alongside them.

  SCFAs shows the floor of the method: with m = 4, all of them above 0.17, SPLINE
  claims 1.06 nulls. DECREASESLOPE returns 4 - i.e. it correctly declines to find
  any signal, and its q-values equal BH's. m = 4 is too few to estimate m0 at
  all, and the right estimator says so rather than guessing.

  K1ROI is the instructive edge case. Thirteen of fourteen p-values are BELOW
  0.05, so the benchmark is low (1.1) and a small m0 is genuinely plausible here
  - but SPLINE returns exactly 0, which asserts no nulls exist and cannot be
  right. DECREASESLOPE returns 14, the opposite extreme: here it is the
  conservative one, and its q-values equal BH's. This is the one dataset of the
  five where the lab default costs power rather than saving it, and it is the
  shape of data - almost everything significant - where that trade is cheap.

  -------------------------------------------------------------------------
  by: Lukas Van Oudenhove  |  KU Leuven, October 2026
  -------------------------------------------------------------------------
  proc_multtest_fdr.sas   v1.3   last modified: 2026/10/01
*****************************************************************************/

/* m=7, 2 p-values > 0.05. benchmark m0 ~ 2.1; DECREASESLOPE m0 = 2.00, exactly on the
   benchmark. (The SAS default spline would give 1.24 - also passing, but lower.) */
data CytokinesT1;
input Raw_P;
datalines;
0.0063
0.0046
0.1097
0.0017
0.0190
0.0025
0.8031
;

ods graphics on;
proc multtest inpvalues=cytokinesT1 fdr pfdr afdr m0=decreaseslope plots=all;
run;

/* m=7, 6 p-values > 0.05. benchmark m0 ~ 6.3; DECREASESLOPE m0 = 3.00. Better than the
   spline's 0.18, which is absurd, but still half the benchmark - this one does NOT fully
   pass. Report plain FDR alongside, or in place of, the q-values below. */
data CytokinesT2;
input Raw_P;
datalines;
0.3913
0.0124
0.2928
0.1349
0.2515
0.0839
0.7543
;
ods graphics on;
proc multtest inpvalues=cytokinesT2 fdr pfdr afdr m0=decreaseslope plots=all;
run;

/* m=14, 12 p-values > 0.05. benchmark m0 ~ 12.6; DECREASESLOPE m0 = 11.00 -> passes.
   (Spline 7.50: passes too, but buys significance the data does not support.) */
data VTROI;
input Raw_P;
datalines;
0.0030
0.0104
0.9732
0.4175
0.0604
0.8024
0.0749
0.1635
0.5736
0.3671
0.5553
0.3952
0.0810
0.8423
;
ods graphics on;
proc multtest inpvalues=VTROI fdr pfdr afdr m0=decreaseslope plots=all;
run;

/* m=14, only 1 p-value > 0.05 so the benchmark is low (~1.1) and a small m0 IS plausible.
   The spline returns m0 = 0, which asserts no nulls exist and cannot be right. DECREASESLOPE
   returns 14.00, so here the lab default is the CONSERVATIVE choice and its q-values equal
   BH's. The one dataset of the five where the new default costs power - cheaply, since
   thirteen of fourteen are significant either way. */
data K1ROI;
input Raw_P;
datalines;
0.0008
0.0013
0.0215
0.0246
0.0003
0.0005
0.0162
0.0422
0.0343
0.0105
0.0026
0.0010
0.0350
0.0699
;
ods graphics on;
proc multtest inpvalues=K1ROI fdr pfdr afdr m0=decreaseslope plots=all;
run;

/* m=4, all 4 p-values > 0.05. benchmark m0 ~ 4.2; DECREASESLOPE m0 = 4.00, i.e. it declines
   to find signal and reproduces plain FDR exactly. (Spline: 1.06, which FAILS.) m=4 is too
   few to estimate m0 at all, and the right estimator says so rather than guessing. */
data SCFAs;
input Raw_P;
datalines;
0.6683
0.8466
0.1787
0.2235
;
ods graphics on;
proc multtest inpvalues=SCFAs fdr pfdr afdr m0=decreaseslope plots=all;
run;