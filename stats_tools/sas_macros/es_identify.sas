/*===========================================================================*
 |  es_identify.sas
 |
 |  Forensic helper: given an effect size that does not match its F and
 |  degrees of freedom, work out which statistic it actually is.
 |
 |  Written for auditing existing results tables -- including tables
 |  produced by the older effect-size script, in which the quantity
 |  reported as "partial eta-squared" was in fact eta-squared.
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
 |  audit -- not only ones using the local adaptation.
 |
 |  ---------------------------------------------------------------------
 |  WHEN THE OLD FORMULA IS ACCIDENTALLY RIGHT
 |  ---------------------------------------------------------------------
 |  The old partial_eta_2 reduces to  qF / (qF + _freq_ - q).  The correct
 |  value is  qF / (qF + DenDF).  They agree IF AND ONLY IF
 |
 |      _freq_ - q  =  DenDF
 |
 |  which is a property of the DESIGN, not of the code:
 |
 |    * NO between-subject factor -- a within-subject time factor crossed
 |      with a covariate, say.  DenDF is observation-driven, so _freq_ - q
 |      lands within a few units of it and the error disappears into the
 |      third decimal.  Such tables audit clean.  This is why the bug
 |      survived: whole projects can validate perfectly.
 |    * A BETWEEN-subject factor in a repeated design.  DenDF is
 |      participant-driven while _freq_ counts observations; they diverge
 |      by the number of repeats per participant.
 |
 |  The published paper's own Case 1 (fev1, 72 patients x 8 hours = 576
 |  observations) demonstrates it: recomputed against the design's real
 |  denominator df, the two WITHIN-subject effects move by 1.2x and the
 |  BETWEEN-subject Drug effect by 6.9x (.029 -> .199).
 |
 |  SO: when auditing an old table, the question is not "was this script
 |  used" but "does the effect being reported cross subjects".  Audit
 |  between-subject rows first.
 |
 |  A caution on eta-squared specifically: the constant-ratio test below
 |  identifies WHICH statistic was reported, but it cannot tell you
 |  whether the eta-squared itself was computed from the right variable.
 |  A copy-paste error that feeds another outcome's variance into
 |  ss_total changes only the constant, so every row stays mutually
 |  consistent and the test passes.  Check the source code too.
 |
 |  ---------------------------------------------------------------------
 |  THE THREE CANDIDATES
 |  ---------------------------------------------------------------------
 |  For an effect with numerator df q, denominator df v, F statistic F,
 |  and N observations in the long-format dataset:
 |
 |    1. Correct partial eta-squared      qF / (qF + v)
 |
 |    2. The old script's "partial_eta_2"  qF / (qF + N - q)
 |       Wrong: it substitutes the number of OBSERVATIONS minus the
 |       NUMERATOR df for the model's DENOMINATOR df.  Note that the MSE
 |       cancels out of that expression entirely, which is why the result
 |       depends only on q, F and N.
 |
 |    3. Eta-squared                       qF x R / N
 |       where R = residual variance / total variance.  R is a property of
 |       the MODEL, not of the effect, so it must come out the same for
 |       every effect in one model.  That gives a decisive test.
 |
 |  ---------------------------------------------------------------------
 |  THE DECISIVE TEST
 |  ---------------------------------------------------------------------
 |  If the reported values are eta-squared, then
 |
 |      reported / (q x F)
 |
 |  must be CONSTANT across every effect within a single model, whatever
 |  the sample size happens to be.  This needs no assumptions at all --
 |  not even N.  Run %es_identify over every effect of one model and
 |  compare the RATIO column.
 |
 |  Agreement to within a few per cent (allowing for rounding in the
 |  published values) means the reported values are eta-squared.
 |  Disagreement by a factor of two or more means the value is something
 |  else and must be checked against the original output.
 |
 |  An implied R above 1 is impossible -- residual variance cannot exceed
 |  total variance -- and is conclusive evidence that the number is not
 |  eta-squared either.
 |
 |  ---------------------------------------------------------------------
 |  USAGE
 |  ---------------------------------------------------------------------
 |  One effect at a time:
 |
 |      %es_identify(numdf = 1, f = 40.45, dendf = 161,
 |                   nobs = 664, reported = 0.028);
 |
 |  NOBS= is optional.  Without it the macro still reports the correct
 |  partial eta-squared and the sample-size-free ratio, but cannot
 |  evaluate the old script's formula or the implied R.
 |
 |  A whole model at once, which is how the test is meant to be used:
 |
 |      data audit;
 |          length effect $32;
 |          infile datalines dsd;
 |          input effect $ numdf dendf f reported;
 |          datalines;
 |      group,1,161,40.45,0.028
 |      time,3,161,66.22,0.130
 |      group*time,3,161,2.44,0.005
 |      ;
 |      run;
 |
 |      %es_identify_ds(data = audit, nobs = 664);
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
%macro es_identify(numdf = , f = , dendf = , nobs = , reported = );

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

        %if %superq(nobs) ne %then %do;
            nobs          = &nobs;
            old_buggy_pe2 = (numdf*f)/(numdf*f + nobs - numdf);
            implied_R     = reported*nobs/(numdf*f);
            put "  Old script's partial_eta_2 (buggy)      : " old_buggy_pe2 8.4;
            put "  Residual/total variance R that the";
            put "    reported value implies, read as eta2  : " implied_R 8.4;
            if implied_R > 1 then
                put "    -> R > 1 is IMPOSSIBLE. Not eta-squared either; check the source output.";
            else if implied_R <= 0 then
                put "    -> R <= 0 is impossible; check the source output.";
            else
                put "    -> plausible; confirm it is the same for every effect in this model.";
        %end;

        put "  Sample-size-free diagnostic";
        put "    reported / (NumDF x F)                : " ratio e10.;
        put "    (must be CONSTANT across all effects in one model";
        put "     if the reported values are eta-squared)";
        put "-------------------------------------------------------------";
    run;

%mend es_identify;


/*---------------------------------------------------------------------------*
 |  %es_identify_ds -- a whole model at once, as a dataset and a table
 |
 |  data=      dataset with one row per effect and the variables
 |             EFFECT (optional), NUMDF, DENDF, F, REPORTED
 |  nobs=      number of observations in the model's long-format dataset;
 |             optional, but needed for the implied-R column
 |  out=       output dataset, default ES_AUDIT
 |  tol=       tolerance for the consistency verdict, as a max/min ratio
 |             of the RATIO column.  Default 1.15, which allows for the
 |             rounding in a published table.
 *---------------------------------------------------------------------------*/
%macro es_identify_ds(data = , nobs = , out = es_audit, tol = 1.15);

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

        %if %superq(nobs) ne %then %do;
            nobs          = &nobs;
            old_buggy_pe2 = (numdf*f)/(numdf*f + nobs - numdf);
            implied_R     = reported*nobs/(numdf*f);
            length R_verdict $24;
            if implied_R > 1        then R_verdict = 'IMPOSSIBLE (R > 1)';
            else if implied_R <= 0  then R_verdict = 'IMPOSSIBLE (R <= 0)';
            else                         R_verdict = 'plausible';
        %end;

        label
            correct_pe2   = 'Correct partial eta-squared'
            reported      = 'Reported value'
            ratio         = 'reported / (NumDF x F)'
            discrepancy   = 'Correct / reported'
        %if %superq(nobs) ne %then %do;
            old_buggy_pe2 = "Old script's partial_eta_2 (buggy)"
            implied_R     = 'Implied residual/total variance R'
            R_verdict     = 'Is that R possible?'
        %end;
        ;
    run;

    /* Consistency of the RATIO column is the decisive test */
    %local _mx _mn _spread _verdict;
    %let _mx = .;  %let _mn = .;  %let _spread = n/a;
    proc sql noprint;
        select max(ratio), min(ratio) into :_mx trimmed, :_mn trimmed from &out;
    quit;

    %if %superq(_mn) = or %superq(_mn) = . %then
        %let _verdict = cannot evaluate (no usable rows in &data);
    %else %if %sysevalf(&_mn > 0) %then %do;
        %let _spread = %sysevalf(&_mx / &_mn);
        %if %sysevalf(&_spread <= &tol) %then
            %let _verdict = CONSISTENT -- the reported values are eta-squared;
        %else
            %let _verdict = INCONSISTENT -- at least one value is not eta-squared;
    %end;
    %else %let _verdict = cannot evaluate (a ratio is zero or negative);

    title  "Effect-size audit";
    title2 "reported / (NumDF x F) spread = &_spread  (tolerance &tol)";
    title3 "Verdict: &_verdict";
    proc print data = &out noobs label;
        format correct_pe2 reported 8.4 ratio e10. discrepancy 8.2
        %if %superq(nobs) ne %then %str(old_buggy_pe2 8.4 implied_R 8.4);
        ;
    run;
    title;

    %put NOTE: (es_identify_ds) ratio spread = &_spread; verdict: &_verdict;
    %put NOTE- Audit written to &out..;

%mend es_identify_ds;


/*===========================================================================*
 |  WORKED EXAMPLE -- the audit of the published MAST models
 |
 |  All three effects come from one model, so if the reported values are
 |  eta-squared then reported/(NumDF x F) must agree across all three.
 |  It does: 6.92e-4, 6.54e-4, 6.83e-4, a spread of 1.06.  They are
 |  eta-squared, reported under the name partial eta-squared.
 *---------------------------------------------------------------------------*
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
 *===========================================================================*/
