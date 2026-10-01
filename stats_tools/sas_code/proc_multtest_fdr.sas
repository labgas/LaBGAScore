/*****************************************************************************
  proc_multtest_fdr.sas

  Multiple-comparison correction with PROC MULTTEST, as a worked example for the
  lab. Each block below feeds one set of raw p-values to

      proc multtest inpvalues=<data> fdr pfdr afdr plots=all;

  which requests three corrections side by side, plus the diagnostic graphics:

      FDR    Benjamini-Hochberg. Assumes every hypothesis could be null
             (pi0 = 1). Always valid, never the most powerful.
      PFDR   Storey's positive FDR. Scales BH by an ESTIMATED number of true
             nulls (m0). More powerful when that estimate is good, and
             anti-conservative when it is not.
      AFDR   Adaptive FDR. Also estimates m0, but by a different route
             (LOWESTSLOPE, after Schweder & Spjotvoll 1982), which is built for
             small numbers of tests and does not need a well-populated upper
             tail of the p-value distribution.
      PLOTS=ALL  the diagnostic plots, including the one that matters here: the
             fitted estimate of the number of true nulls.

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

  1. Try a different m0 estimator. The option is NTRUENULL=, with M0= as a
     documented alias; PTRUENULL= takes a PROPORTION instead of a count, and
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
     conclude that two adjustments disagree about the data:

         PFDR            SPLINE first; if the estimate is nonpositive, or if the
                         slope of the spline at the last lambda exceeds 0.1 times
                         the range of the fitted spline values, BOOTSTRAP instead
         ADAPTIVEFDR     LOWESTSLOPE
         ADAPTIVEHOLM    DECREASESLOPE
         ADAPTIVEHOCHBERG  DECREASESLOPE

     LOWESTSLOPE and DECREASESLOPE are the ones to reach for at small m: they
     read m0 off the slope of the ordered p-values and need no well-populated
     upper tail, which is exactly what SPLINE and BOOTSTRAP lack when only a
     handful of p-values exceed 0.05.

  2. Or report ADAPTIVEFDR (alias AFDR) instead - but only if ITS m0 passes the
     same check. It defaults to LOWESTSLOPE for that reason, and in the five
     examples below it passes every time PFDR fails.

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

  Computed with LaBGAScore_Storey_FDR, which reproduces PROC MULTTEST's PFDR
  spline exactly (verified on CytokinesT1/T2 to five decimals):

    dataset        m   #p>.05  benchmark m0   PFDR m0   AFDR m0   verdict
    CytokinesT1    7      2        2.1          1.24      6.00    both pass
    CytokinesT2    7      6        6.3          0.18      5.00    PFDR FAILS
    VTROI         14     12       12.6          7.50     13.00    both pass
    K1ROI         14      1        1.1          0.00     14.00    PFDR m0 = 0
    SCFAs          4      4        4.2          1.06      4.00    PFDR FAILS

  CytokinesT2 is the clearest failure: six of seven p-values exceed 0.05, yet
  PFDR puts the number of true nulls at 0.18. AFDR says 5, which is credible.

  SCFAs shows the floor of the method: with m = 4, all of them above 0.17, PFDR
  claims 1.06 nulls. m = 4 is too few to estimate m0 at all - report FDR.

  K1ROI is the instructive edge case. Thirteen of fourteen p-values are BELOW
  0.05, so the benchmark is low (1.1) and a small m0 is genuinely plausible here
  - but PFDR returns exactly 0, which asserts no nulls exist and cannot be right.
  AFDR's 14 is the opposite extreme, over-conservative given the data. When the
  two bracket the answer this widely, say which you used and why.

  -------------------------------------------------------------------------
  by: Lukas Van Oudenhove  |  KU Leuven, October 2026
  -------------------------------------------------------------------------
  proc_multtest_fdr.sas   v1.1   last modified: 2026/10/01
*****************************************************************************/

/* m=7, 2 p-values > 0.05. benchmark m0 ~ 2.1; PFDR 1.24, AFDR 6.00 -> both pass */
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
proc multtest inpvalues=cytokinesT1 fdr pfdr afdr plots=all;
run;

/* m=7, 6 p-values > 0.05. benchmark m0 ~ 6.3; PFDR 0.18 FAILS the check, AFDR 5.00 passes -> report AFDR */
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
proc multtest inpvalues=cytokinesT2 fdr pfdr afdr plots=all;
run;

/* m=14, 12 p-values > 0.05. benchmark m0 ~ 12.6; PFDR 7.50, AFDR 13.00 -> both pass, PFDR notably lower */
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
proc multtest inpvalues=VTROI fdr pfdr afdr plots=all;
run;

/* m=14, only 1 p-value > 0.05 so the benchmark is low (~1.1) and a small m0 IS plausible -
   but PFDR returns m0 = 0, which asserts no nulls exist and cannot be right. AFDR 14.00 is the
   other extreme. State which you used. */
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
proc multtest inpvalues=K1ROI fdr pfdr afdr plots=all;
run;

/* m=4, all 4 p-values > 0.05. benchmark m0 ~ 4.2; PFDR 1.06 FAILS, AFDR 4.00 = FDR.
   m=4 is too few to estimate m0 at all -> report plain FDR */
data SCFAs;
input Raw_P;
datalines;
0.6683
0.8466
0.1787
0.2235
;
ods graphics on;
proc multtest inpvalues=SCFAs fdr pfdr afdr plots=all;
run;