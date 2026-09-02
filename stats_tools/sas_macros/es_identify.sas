/*===========================================================================*
 |  es_identify.sas
 |
 |  Forensic helper: given an effect size that does not match its F and
 |  degrees of freedom, work out which statistic it actually is and how it
 |  was computed.
 |
 |  Companion file:  mixed_effectsize.sas  -- computes effect sizes
 |  correctly for new analyses.
 |
 |  ---------------------------------------------------------------------
 |  PROVENANCE -- WHAT YOU ARE ACTUALLY AUDITING
 |  ---------------------------------------------------------------------
 |  The old script derives from a published SAS Users Group paper:
 |
 |      Tippey KG, Longnecker MT.  An Ad Hoc Method for Computing
 |      Pseudo-Effect Size for Mixed Models.  Texas A&M University.
 |
 |  The line  ss_error = mse * (_freq_ - numdf);  is verbatim in that
 |  paper's macro appendix.  So this is a defect in the method as
 |  published, and tables from ANY project built on it are candidates for
 |  audit -- not only ones using a local adaptation.
 |
 |  ---------------------------------------------------------------------
 |  WHAT WENT WRONG IN PRACTICE
 |  ---------------------------------------------------------------------
 |  The old script writes three statistics into one dataset -- eta_2,
 |  omega_2 and partial_eta_2 -- so the first hypothesis is always that
 |  the wrong COLUMN was reported.  In the ME/CFS case-control analysis
 |  that turned out NOT to be what happened.  The reported quantity was
 |  partial_eta_2 itself: the name was right and the arithmetic was
 |  wrong, which is the harder failure to spot.
 |
 |  Two things compound in that arithmetic:
 |
 |    1. _freq_ - numdf stands in for the model's DENOMINATOR df.
 |    2. _freq_ is the number of rows PROC MEANS saw in the OUTPM=
 |       dataset -- that is, observations READ, not observations USED.
 |       An OUTPM= dataset keeps every row the procedure read, including
 |       those dropped for a missing response or covariate.
 |
 |  In that analysis PROC MIXED read 1384 rows for every model and used
 |  between 483 and 1312 of them, so _freq_ was 1384 throughout while the
 |  denominator df ranged from 124 to 161.  Published partial eta-squared
 |  values were between 3x and 9x too small.
 |
 |  ---------------------------------------------------------------------
 |  THE DECISIVE TEST -- SOLVE FOR THE IMPLIED _FREQ_
 |  ---------------------------------------------------------------------
 |  If a reported value r came from the buggy formula
 |
 |      r = qF / (qF + f - q)
 |
 |  then inverting it recovers the f that produced it:
 |
 |      implied_freq = qF(1 - r)/r + q
 |
 |  Compute this for every effect in the table.  If the values cluster on
 |  ONE number, and that number matches the row count of the source
 |  dataset (the OUTPM= dataset, or "Number of Observations Read" in the
 |  PROC MIXED output), the reported values are the buggy partial_eta_2
 |  and the diagnosis is closed.  The clustering is itself the evidence:
 |  f is a property of the dataset, so it cannot vary between effects.
 |
 |  Worked example -- the ME/CFS models.  Twenty published values across
 |  seven models implied f between 1313 and 1485, median 1381, against a
 |  printed Number of Observations Read of 1384 in every one.  The spread
 |  is entirely two-significant-figure rounding in the published values.
 |
 |  ---------------------------------------------------------------------
 |  THE SECONDARY TEST -- AND WHY IT IS NOT DECISIVE ON ITS OWN
 |  ---------------------------------------------------------------------
 |  If the reported values were eta_2 instead, then
 |
 |      r / (q x F)
 |
 |  would be constant across every effect within one model, because it
 |  equals (residual/total variance) / N, both properties of the model
 |  rather than of the effect.
 |
 |  *** THIS TEST CANNOT SEPARATE eta_2 FROM THE BUGGY partial_eta_2. ***
 |
 |  When f is much larger than qF -- which is the normal case, since f
 |  counts observations -- the buggy formula behaves like
 |
 |      qF / (qF + f - q)  ~=  qF / f
 |
 |  which is also proportional to qF.  So BOTH candidates produce a
 |  near-constant ratio and the test passes either way.  In the ME/CFS
 |  audit it passed at a spread of 1.14 against a tolerance of 1.15, and
 |  the wrong conclusion survived several rounds of review because of it.
 |
 |  Use the ratio test only to rule out a value that fits NEITHER
 |  candidate.  Use the implied-_freq_ test to decide between them.
 |
 |  A further limitation of both tests: neither can detect an eta_2
 |  computed from the WRONG VARIABLE.  Feeding another outcome's variance
 |  into ss_total changes only the shared constant, so every row stays
 |  mutually consistent and both tests pass.  Read the source code as
 |  well; that is how the RMSSD/HR variance swap in the ME/CFS script was
 |  found.
 |
 |  ---------------------------------------------------------------------
 |  WHEN THE OLD FORMULA IS ACCIDENTALLY RIGHT
 |  ---------------------------------------------------------------------
 |  The old partial_eta_2 reduces to  qF / (qF + f - q).  The correct
 |  value is  qF / (qF + DenDF).  They agree IF AND ONLY IF
 |
 |      f - q  =  DenDF
 |
 |  which needs TWO independent things to hold: f must count only the
 |  rows the model FITTED, and DenDF must be the observation-scale
 |  residual df.
 |
 |  The second is decided by the DF METHOD in force, NOT by whether the
 |  effect crosses subjects:
 |
 |    CONTAINMENT, falling through -- a RANDOM statement is present and
 |      no random effect contains the fixed effect
 |        -> N - rank(XZ) for EVERY effect.  AGREES.  Whole projects can
 |           validate perfectly, which is why the bug survived.
 |    BETWEEN-WITHIN with a residual variance term (CS, VC, AR(1))
 |        -> stratifies.  Agrees for within-subject effects, FAILS for
 |           between-subject ones.
 |    BETWEEN-WITHIN under type = UN
 |        -> participant scale for EVERY effect.  FAILS for all.
 |    KR / KR2 / SATTERTH
 |        -> fractional.  FAILS.
 |
 |  The published paper's own Case 1 (fev1, 72 patients x 8 hours = 576
 |  observations) is the SECOND case: recomputed against the design's
 |  real denominator df, the two WITHIN-subject effects move by 1.2x and
 |  the BETWEEN-subject Drug effect by 6.9x (.029 -> .199).  Do not
 |  generalise from it.
 |
 |  SO: when auditing an old table, the question is not "was this script
 |  used", nor "does the effect cross subjects", but "which DF METHOD was
 |  in force, and did f count only fitted rows".
 |
 |  ---------------------------------------------------------------------
 |  USAGE
 |  ---------------------------------------------------------------------
 |  One effect at a time:
 |
 |      %es_identify(numdf = 1, f = 40.45, dendf = 161,
 |                   nread = 1384, reported = 0.028);
 |
 |  NREAD= is the row count of the dataset PROC MEANS ran on -- normally
 |  "Number of Observations Read" from the PROC MIXED output, NOT
 |  "Number of Observations Used".  Optional, but it is what makes the
 |  decisive test possible.  NOBS= (observations used) is optional too
 |  and only feeds the implied residual/total variance ratio.
 |
 |  A whole model at once, which is how the test is meant to be used:
 |
 |      data audit;
 |          length effect $32;
 |          infile datalines dsd;
 |          input effect $ numdf dendf f reported;
 |          datalines;
 |      Time_VAS_U,3,161,66.22,0.13
 |      Group,1,161,40.45,0.028
 |      Time_VAS_U*Group,3,161,2.44,0.005
 |      ;
 |      run;
 |
 |      %es_identify_ds(data = audit, nread = 1384, nobs = 656);
 |
 |  ---------------------------------------------------------------------
 |  Written for the Laboratory for Brain-Gut Axis Studies (LaBGAS),
 |  KU Leuven.
 |
 |  NOT YET RUN IN SAS -- the arithmetic has been verified independently,
 |  but please confirm the output against a known case before relying on it.
 *===========================================================================*/


/*---------------------------------------------------------------------------*
 |  %es_identify -- one effect, printed to the log
 *---------------------------------------------------------------------------*/
%macro es_identify(numdf = , f = , dendf = , nread = , nobs = , reported = );

    %if %superq(numdf) = or %superq(f) = or %superq(dendf) = or %superq(reported) = %then %do;
        %put ERROR: (es_identify) NUMDF=, F=, DENDF= and REPORTED= are all required.;
        %return;
    %end;

    data _null_;
        numdf = &numdf;  f = &f;  dendf = &dendf;  reported = &reported;

        correct_pe2 = (numdf*f)/(numdf*f + dendf);
        ratio       = reported/(numdf*f);

        put "-------------------------------------------------------------";
        put "  F(" numdf best. ", " dendf best. ") = " f best. "    reported = " reported best.;
        put "-------------------------------------------------------------";
        put "  Correct partial eta-squared             : " correct_pe2 8.4;
        put "  Correct / reported                      : " (correct_pe2/reported) 8.2 "x";

        /* --- decisive test: what _freq_ would the buggy formula have needed? --- */
        if 0 < reported < 1 then do;
            implied_freq = numdf*f*(1 - reported)/reported + numdf;
            put " ";
            put "  DECISIVE TEST -- implied _freq_";
            put "    qF(1-r)/r + q                        : " implied_freq 10.0;
            %if %superq(nread) ne %then %do;
                nread = &nread;
                put "    rows actually read by the procedure   : " nread 10.0;
                put "    ratio                                 : " (implied_freq/nread) 8.3;
                if 0.90 <= implied_freq/nread <= 1.10 then
                    put "    -> MATCH. The reported value is the buggy partial_eta_2.";
                else
                    put "    -> no match; try the ratio test below across all effects.";
                buggy_pe2 = (numdf*f)/(numdf*f + nread - numdf);
                put "    buggy partial_eta_2 recomputed        : " buggy_pe2 8.4;
            %end;
            %else %do;
                put "    (supply NREAD= -- Number of Observations READ -- to complete this test)";
            %end;
        end;

        /* --- secondary test: eta-squared --- */
        put " ";
        put "  SECONDARY TEST -- eta-squared";
        put "    reported / (NumDF x F)                : " ratio e10.;
        put "    (must be CONSTANT across all effects of one model if the";
        put "     reported values are eta-squared -- but see the header: the";
        put "     buggy partial_eta_2 also gives a near-constant ratio, so";
        put "     this CANNOT separate the two on its own)";
        %if %superq(nobs) ne %then %do;
            nobs      = &nobs;
            implied_R = reported*nobs/(numdf*f);
            put "    implied residual/total variance R     : " implied_R 8.4;
            if implied_R > 1 then
                put "    -> R > 1 is IMPOSSIBLE. Not eta-squared.";
            else if implied_R <= 0 then
                put "    -> R <= 0 is impossible; check the source output.";
        %end;
        put "-------------------------------------------------------------";
    run;

%mend es_identify;


/*---------------------------------------------------------------------------*
 |  %es_identify_ds -- a whole model at once, as a dataset and a table
 |
 |  data=    dataset with one row per effect and the variables
 |           EFFECT (optional), NUMDF, DENDF, F, REPORTED
 |  nread=   rows in the dataset PROC MEANS ran on = "Number of
 |           Observations Read". Optional but strongly recommended: it is
 |           what makes the decisive test possible.
 |  nobs=    observations USED. Optional; feeds the implied-R column only.
 |  out=     output dataset, default ES_AUDIT
 |  tol=     tolerance for the eta-squared consistency verdict, as a
 |           max/min ratio of the RATIO column. Default 1.15.
 *---------------------------------------------------------------------------*/
%macro es_identify_ds(data = , nread = , nobs = , out = es_audit, tol = 1.15);

    %if %superq(data) = %then %do;
        %put ERROR: (es_identify_ds) DATA= is required.;
        %return;
    %end;
    %else %if not %sysfunc(exist(&data)) %then %do;
        %put ERROR: (es_identify_ds) Dataset &data does not exist.;
        %return;
    %end;

    data &out;
        set &data;

        correct_pe2 = (numdf*f)/(numdf*f + dendf);
        ratio       = reported/(numdf*f);
        discrepancy = correct_pe2 / reported;

        if 0 < reported < 1 then implied_freq = numdf*f*(1 - reported)/reported + numdf;

        %if %superq(nread) ne %then %do;
            nread     = &nread;
            buggy_pe2 = (numdf*f)/(numdf*f + nread - numdf);
            freq_ratio = implied_freq / nread;
        %end;
        %if %superq(nobs) ne %then %do;
            nobs      = &nobs;
            implied_R = reported*nobs/(numdf*f);
            length R_verdict $24;
            if implied_R > 1        then R_verdict = 'IMPOSSIBLE (R > 1)';
            else if implied_R <= 0  then R_verdict = 'IMPOSSIBLE (R <= 0)';
            else                         R_verdict = 'plausible';
        %end;

        label
            correct_pe2  = 'Correct partial eta-squared'
            reported     = 'Reported value'
            implied_freq = 'Implied _freq_  qF(1-r)/r + q'
            ratio        = 'reported / (NumDF x F)'
            discrepancy  = 'Correct / reported'
        %if %superq(nread) ne %then %do;
            buggy_pe2    = "Buggy partial_eta_2 recomputed at _freq_ = &nread"
            freq_ratio   = 'Implied _freq_ / rows read'
        %end;
        %if %superq(nobs) ne %then %do;
            implied_R    = 'Implied residual/total variance R'
            R_verdict    = 'Is that R possible?'
        %end;
        ;
    run;

    /*-- Decisive test: do the implied _freq_ values cluster on NREAD= ? --*/
    %local _fmin _fmax _fmed _verdict _spread _rmx _rmn;
    %let _verdict = ;
    %if %superq(nread) ne %then %do;
        proc sql noprint;
            select min(freq_ratio), max(freq_ratio), median(implied_freq)
              into :_fmin trimmed, :_fmax trimmed, :_fmed trimmed
              from &out where implied_freq is not missing;
        quit;
        %if %superq(_fmin) ne and %superq(_fmin) ne . %then %do;
            %if %sysevalf(&_fmin >= 0.85) and %sysevalf(&_fmax <= 1.15) %then
                %let _verdict = BUGGY PARTIAL_ETA_2 -- implied _freq_ (median &_fmed) matches the &nread rows read;
        %end;
    %end;

    /*-- Secondary test: is the eta-squared ratio constant? --*/
    %let _spread = n/a;
    proc sql noprint;
        select max(ratio), min(ratio) into :_rmx trimmed, :_rmn trimmed from &out;
    quit;
    %if %superq(_rmn) ne and %superq(_rmn) ne . %then %do;
        %if %sysevalf(&_rmn > 0) %then %let _spread = %sysevalf(&_rmx / &_rmn);
    %end;

    %if %superq(_verdict) = %then %do;
        %if %superq(_spread) ne n/a and %sysevalf(&_spread <= &tol) %then
            %let _verdict = consistent with eta-squared (ratio spread &_spread) -- but supply NREAD= to rule out the buggy partial_eta_2;
        %else
            %let _verdict = INCONSISTENT -- fits neither candidate; check the source output;
    %end;

    title  "Effect-size audit";
    title2 "Decisive test: implied _freq_ = qF(1-r)/r + q, compared with the rows READ";
    title3 "Verdict: &_verdict";
    title4 "Secondary: reported/(NumDF x F) spread = &_spread (a near-constant ratio does NOT distinguish eta-squared from the buggy partial_eta_2)";
    proc print data = &out noobs label;
        format correct_pe2 reported 8.4 implied_freq 10.0 ratio e10. discrepancy 8.2
        %if %superq(nread) ne %then %str(buggy_pe2 8.4 freq_ratio 8.3);
        %if %superq(nobs)  ne %then %str(implied_R 8.4);
        ;
    run;
    title;

    %put NOTE: (es_identify_ds) verdict: &_verdict;
    %put NOTE- Audit written to &out..;

%mend es_identify_ds;


/*===========================================================================*
 |  WORKED EXAMPLE -- the audit of the published MAST self-reported stress
 |  model, using the true F values from the SAS output.
 |
 |  Implied _freq_ comes out at 1405 / 1332 / 1460 against a printed
 |  Number of Observations Read of 1384.  Verdict: the published values
 |  are the buggy partial_eta_2, not eta-squared.
 |
 |  Note that the eta-squared ratio test ALSO passes here, at a spread of
 |  1.14 -- which is precisely why it must not be used on its own.
 *---------------------------------------------------------------------------*
data mast_stress;
    length effect $32;
    infile datalines dsd;
    input effect $ numdf dendf f reported;
    datalines;
Time_VAS_U,3,161,66.22,0.13
Group,1,161,40.45,0.028
Time_VAS_U*Group,3,161,2.44,0.005
;
run;

%es_identify_ds(data = mast_stress, nread = 1384, nobs = 656);
 *===========================================================================*/
