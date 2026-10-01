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

  1. Try a different m0 estimator for PFDR. SAS's PFDR default is the SPLINE
     method, falling back to BOOTSTRAP; DECREASESLOPE and LOWESTSLOPE are the
     alternatives, and are steadier when the number of tests is small. Set it
     with the NTRUENULL= option on the PROC MULTTEST statement. CHECK THE
     SYNTAX against the SAS documentation for your installed version: nothing in
     this folder has been executed in SAS (see sas_macros/README.md).

  2. Or report AFDR instead - but only if ITS m0 passes the same check. AFDR's
     estimator behaves far better at small m, and in the five examples below it
     passes every time PFDR fails.

  3. Or report plain FDR (Benjamini-Hochberg). Always defensible, and the honest
     answer when m0 simply is not estimable - which it is not when only a handful
     of p-values sit above 0.05.

  Do NOT respond to a failed check by specifying m0 by hand unless you have
  external grounds for the number. Specifying m0 = m is exactly BH; anything
  lower buys significance by assumption rather than from the data.

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