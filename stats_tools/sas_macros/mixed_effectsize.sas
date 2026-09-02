/*===========================================================================*
 |  mixed_effectsize.sas
 |
 |  Effect sizes for fixed effects in models fitted with PROC MIXED --
 |  MARGINAL models (REPEATED statement, population-average) and true
 |  MIXED models (RANDOM statement, subject-specific) alike.
 |
 |  THE FULL DOCUMENTATION IS IN  sas_macros/README.md.  This header is
 |  the working reference -- usage, parameters, output, and the choices
 |  you have to make before calling the macro.  The README carries the
 |  rest: why this macro exists and what defect in the published method
 |  it replaces, the observations-read-versus-used problem, the
 |  cross-validation against R and Python, the correspondence with
 |  PROC GLM, the caveats to state in a Methods section, and the
 |  references.
 |
 |  Companion file:  es_identify.sas -- forensic helper for auditing
 |  effect sizes in tables produced by other code.
 |
 |  ---------------------------------------------------------------------
 |  USAGE
 |  ---------------------------------------------------------------------
 |  Marginal model (REPEATED):
 |
 |      proc mixed data = mydata;
 |          class id group time;
 |          model y = group|time / outpm = outpm_y ddfm = kr2;
 |          repeated time / subject = id type = un;
 |          ods output tests3 = tests3_y;
 |      run;
 |
 |      %mixed_effectsize(tests3 = tests3_y,
 |                        model  = MARGINAL,
 |                        resids = outpm_y,
 |                        dv     = y,
 |                        out    = es_y,
 |                        label  = %str(Self-reported stress, MAST));
 |
 |  True mixed model (RANDOM), marginal eta-squared -- the default:
 |
 |      proc mixed data = mydata;
 |          class id group time;
 |          model y = group|time / outpm = outpm_y outp = outp_y ddfm = kr2;
 |          random intercept / subject = id;
 |          ods output tests3 = tests3_y;
 |      run;
 |
 |      %mixed_effectsize(tests3 = tests3_y,
 |                        model  = MIXED,
 |                        resids = outpm_y,        <-- OUTPM, not OUTP
 |                        dv     = y);
 |
 |  True mixed model, within-subject (conditional) eta-squared:
 |
 |      %mixed_effectsize(tests3    = tests3_y,
 |                        model     = MIXED,
 |                        residtype = CONDITIONAL,
 |                        resids    = outp_y,      <-- OUTP
 |                        dv        = y);
 |
 |  With the degrees-of-freedom diagnostics switched on:
 |
 |      proc mixed data = mydata;
 |          ...
 |          ods output tests3 = t3  modelinfo = mi  dimensions = dm;
 |      run;
 |
 |      %mixed_effectsize(tests3 = t3, resids = outpm_y, dv = y,
 |                        modelinfo = mi, dimensions = dm,
 |                        nsubjects = 164, between = Group Centered_BMI);
 |
 |  Partial eta-squared and its CI need only TESTS3.  RESIDS= and DV= are
 |  needed only for eta-squared.
 |
 |  ---------------------------------------------------------------------
 |  PARAMETERS
 |  ---------------------------------------------------------------------
 |  tests3     Dataset from ODS OUTPUT TESTS3=.  REQUIRED.  Must contain
 |             Effect, NumDF, DenDF, FValue (ProbF optional).
 |  model      MARGINAL (default) for a REPEATED-only model, MIXED for a
 |             model with a RANDOM statement.  Controls validation and the
 |             wording of the output; the partial eta-squared arithmetic
 |             is identical either way.
 |  residtype  MARGINAL (default) or CONDITIONAL.  Says which residuals
 |             the supplied dataset holds, and hence how eta-squared is
 |             named and interpreted.  Refused when MODEL = MARGINAL,
 |             where OUTP= and OUTPM= are the same thing.
 |  resids     Residuals dataset: OUTPM= when RESIDTYPE = MARGINAL, OUTP=
 |             when CONDITIONAL.  Optional; omit if only partial
 |             eta-squared is wanted.
 |  outpm      Deprecated alias for RESIDS=, kept so older calls still run.
 |  dv         Dependent variable in the residuals dataset.  Required if
 |             RESIDS= is supplied.
 |  resid      Residual variable.  Default Resid, correct for both OUTPM=
 |             and OUTP=.
 |  out        Output dataset.  Default EFFECTSIZES.
 |  alpha      Two-sided alpha for the partial eta-squared CI.  Default
 |             0.05.
 |  label      Free-text model label carried into the output.  Optional.
 |  print      Y (default) prints a formatted table; N suppresses it.
 |
 |  --- eta-squared method (does NOT affect partial eta-squared) --------
 |  eta2_method  FROM_F (default) reconstructs SS_effect from the F test
 |             as NumDF * FValue * (SS_error / DenDF).  Cheap, but only
 |             exact in a fixed-effects ANOVA -- see CHOICE 3 below.
 |             DIRECT refits the same mean model with PROC GLM and takes
 |             real Type III sums of squares.  Slower, but a genuine
 |             variance decomposition.  Requires DATA=, CLASS= and FIXED=.
 |  data       Input dataset for ETA2_METHOD = DIRECT.  Use the dataset
 |             the PROC MIXED model was fitted to.
 |  class      CLASS variables for the DIRECT refit, space separated.
 |  fixed      Right-hand side of the MODEL statement for the DIRECT
 |             refit -- the fixed effects only, exactly as written in
 |             PROC MIXED.
 |
 |  --- degrees-of-freedom diagnostics (all optional) -------------------
 |  modelinfo  Dataset from ODS OUTPUT MODELINFO=.  Supplies the
 |             covariance structure, the df method actually used, the
 |             residual variance method and the subject effect, so the
 |             macro tailors its df checks instead of guessing.
 |  dimensions Dataset from ODS OUTPUT DIMENSIONS=.  Supplies Columns in
 |             Z (whether a RANDOM statement is active) and the subject
 |             count -- but see NSUBJECTS=.
 |  nsubjects  Number of subjects.  Give it explicitly whenever you want
 |             the between-subject df check.  It is needed because the
 |             Dimensions subject count is not always usable: a RANDOM
 |             effect with no SUBJECT= option spans subjects, so SAS
 |             reports Subjects = 1 and Max Obs per Subject = N.  The
 |             macro detects that and asks for this parameter.
 |  between    Space-separated effect names that are BETWEEN-subject,
 |             exactly as they appear in TESTS3.  With NSUBJECTS= this
 |             turns the df check from a hint into a statement: a
 |             between-subject effect cannot have more denominator df
 |             than there are subjects.
 |
 |  ---------------------------------------------------------------------
 |  OUTPUT
 |  ---------------------------------------------------------------------
 |  One row per effect in TESTS3.  Every variable carries a full label,
 |  so nothing can be reported under the wrong heading.
 |
 |    partial_eta2              (NumDF*F) / (NumDF*F + DenDF)
 |    partial_eta2_lcl/_ucl     CI by inverting the noncentral F
 |    partial_epsilon2          NumDF(F-1) / (NumDF*F + DenDF), less
 |                              positively biased, same inputs
 |    cohen_f2                  partial_eta2 / (1 - partial_eta2)
 |    eta2_approx               with RESIDS= and DV=; named
 |      or eta2_cond_approx     eta2_cond_approx under RESIDTYPE=CONDITIONAL
 |    ss_effect ss_error        with RESIDS=
 |      ss_total mse n_obs
 |      nobs_dendf
 |    n_subjects dendf_ratio    with NSUBJECTS=
 |    is_between                with BETWEEN=
 |    eta2_source               FROM_F | GLM_TYPE3 | UNMATCHED
 |    model_label model_type resid_type
 |
 |  OMEGA-SQUARED IS DELIBERATELY NOT COMPUTED -- LaBGAS does not report
 |  it, in either its partial or its classical form.
 |
 |  WHAT TO REPORT: partial_eta2, because a reader can recompute it from
 |  the F and df printed in the same table, which is what makes a results
 |  table checkable.  Optionally partial_epsilon2 beside it.  Eta-squared
 |  cannot be recovered from F and df, so report it only ALONGSIDE
 |  partial eta-squared, never instead, and name it correctly.  For a
 |  two-group contrast a standardised mean difference with a CI (ESTIMATE
 |  or LSMESTIMATE with a sensible divisor) is more interpretable than any
 |  variance-explained measure.
 |
 |  ---------------------------------------------------------------------
 |  CHOICE 1 -- THE DENOMINATOR-DF METHOD  (the most consequential one)
 |  ---------------------------------------------------------------------
 |  Everything here depends on DenDF, and SAS picks one of two defaults
 |  FOR you unless you name a method: CONTAIN when a RANDOM statement is
 |  present, BETWITHIN for REPEATED-only.  Those two are on different
 |  SCALES, so the same data and the same F can give a partial
 |  eta-squared several times smaller in one model than in the other.
 |
 |     ****************************************************************
 |     *  LaBGAS STANDARD:  ddfm = kr2                                 *
 |     *                                                              *
 |     *     model y = <effects> / ddfm = kr2 outpm = <out>;           *
 |     *                                                              *
 |     *  On the MODEL statement of every mixed or marginal model,     *
 |     *  unless there is a stated reason not to.  Fall back in this   *
 |     *  order, and RECORD WHICH WAS USED -- effect sizes depend on   *
 |     *  it:                                                         *
 |     *     kr2  ->  kr        (SAS/STAT older than 12.1 has no KR2)  *
 |     *          ->  satterth  (KR will not converge, or is far too   *
 |     *                         slow)                                *
 |     ****************************************************************
 |
 |    ddfm = kr2   Kenward-Roger with the Kenward & Roger (2009)
 |                 precision estimator.  THE STANDARD.  Inflates the
 |                 covariance matrix of the fixed effects to allow for
 |                 the covariance parameters having been ESTIMATED, and
 |                 derives a Satterthwaite-type df from it, so it changes
 |                 both F and DenDF.
 |    ddfm = kr    Kenward & Roger (1997).  IDENTICAL to kr2 for a
 |                 covariance structure linear in its parameters -- UN,
 |                 CS, VC, TOEP, the LaBGAS common case.  They differ
 |                 only for NONLINEAR structures (AR(1), ARH(1), CSH,
 |                 TOEPH, SP()).  Use only where kr2 is unavailable.
 |    ddfm = kr(firstorder)   subsumed by kr2; no reason to choose it.
 |    ddfm = satterth   df adjustment only, no covariance correction.
 |                 Acceptable fallback.
 |    ddfm = contain | betwithin | residual   the SAS defaults.  Fine in
 |                 large balanced designs and fine under type=UN, but
 |                 they do not allow for the covariance parameters having
 |                 been estimated, and they leave the choice to SAS.
 |    empirical = mbn | morel   goes on the PROC MIXED statement, not
 |                 MODEL.  Consider it INSTEAD of KR when the covariance
 |                 STRUCTURE itself is doubtful.
 |
 |  Three conditions on KR / KR2:
 |    * defined for METHOD = REML (the PROC MIXED default).  A
 |      likelihood-ratio test on FIXED effects needs ML, where KR does
 |      not apply -- do not mix the two in one table.
 |    * if a variance component is estimated at the boundary (0 in
 |      Covariance Parameter Estimates), the adjustment rests on a
 |      Hessian at that boundary.  Check before trusting the df.
 |    * a RANDOM effect with no SUBJECT= spans subjects, so SAS cannot
 |      block the problem.  Fix the blocking before reaching for KR.
 |
 |  UNDER type = UN, DO NOT "FIX" THE DEFAULT.  BETWITHIN prints the same
 |  DenDF for every effect there and that is very nearly right -- two df
 |  out of 160 against the exact multivariate values, at most 0.006 on
 |  partial eta-squared.  KR2 still earns its place, for the covariance
 |  correction rather than the df.
 |
 |  THE MACRO RUNS FINE ON ANY OF THESE, defaults included.  Nothing about
 |  the df method can stop it -- the df diagnostics are NOTE and WARNING
 |  text only.  What changes is the INTERPRETATION: partial eta-squared is
 |  conditioned on DenDF and has no df-free value, so report the df method
 |  beside the effect size.
 |
 |  Why, in full, and the four regimes in which the old published formula
 |  does and does not fail: README.md, "Choose a denominator-df method".
 |
 |  ---------------------------------------------------------------------
 |  CHOICE 2 -- WHICH RESIDUALS  (only affects eta-squared)
 |  ---------------------------------------------------------------------
 |  PARTIAL eta-squared is the same either way: it depends only on NumDF,
 |  DenDF and F, and DenDF already reflects the covariance structure that
 |  was fitted.  ETA-SQUARED changes, because "the residual" differs:
 |
 |    RESIDTYPE = MARGINAL     y - Xb,     from OUTPM=.  SS_error still
 |      (the default)          contains between-subject variance, so
 |                             eta-squared is a proportion of TOTAL
 |                             observed variance -- the honest analogue of
 |                             classical eta-squared, comparable with
 |                             fixed-effects studies.
 |    RESIDTYPE = CONDITIONAL  y - Xb - Zu, from OUTP=.  The random
 |                             effects have absorbed the between-subject
 |                             variance, so the result is a WITHIN-SUBJECT
 |                             proportion.  Systematically larger, and NOT
 |                             comparable with eta-squared elsewhere.
 |
 |  The outputs are named differently (ETA2_APPROX versus
 |  ETA2_COND_APPROX) so the two cannot be confused in a results table.
 |  Ask for CONDITIONAL only when the within-subject proportion is what
 |  you mean to report, say so in the Methods, and report the variance
 |  components beside it.
 |
 |  ---------------------------------------------------------------------
 |  CHOICE 3 -- HOW ETA-SQUARED IS FORMED  (ETA2_METHOD)
 |  ---------------------------------------------------------------------
 |  FROM_F reconstructs SS_effect as NumDF * F * MSE.  That is exact only
 |  in a fixed-effects ANOVA.  Here it can be out by several fold in
 |  EITHER direction, and the error is NOT predictable from the degrees of
 |  freedom -- the macro warns on every row and NOBS_DENDF diagnoses only
 |  half of the problem.  Use DIRECT when eta-squared matters, or report
 |  partial eta-squared and partial epsilon-squared, which are unaffected.
 |  Measured errors and the reasoning: README.md, "Caveats to state in a
 |  Methods section", item 2.
 |
 |  PROC GLM under DIRECT ignores the repeated-measures covariance
 |  structure entirely.  Only its SUMS OF SQUARES are used; its F and p
 |  values are NOT valid for these designs and are never reported.  F, df
 |  and partial eta-squared always come from PROC MIXED.
 |
 |  ---------------------------------------------------------------------
 |  KEY REFERENCES  (full list in the README)
 |  ---------------------------------------------------------------------
 |  Edwards LJ et al.  An R2 statistic for fixed effects in the linear
 |      mixed model.  Stat Med 2008;27(29):6137-6157.  Formalises the
 |      quantity this macro computes as R-squared-beta.
 |  Steiger JH.  Beyond the F test.  Psychol Methods 2004;9(2):164-182.
 |      The noncentral-F confidence interval.
 |  Kenward MG, Roger JH.  Biometrics 1997;53:983-997, and Comput Stat
 |      Data Anal 2009;53:2583-2595.  KR and KR2.
 |  (Verify page numbers before citing.)
 |
 |  ---------------------------------------------------------------------
 |  Written for the Laboratory for Brain-Gut Axis Studies (LaBGAS),
 |  KU Leuven.  Supersedes "Calculation of effect size following marginal
 |  linear mixed models" (B. Dalile, 3 December 2021).
 |
 |  RUN IN SAS: the default FROM_F path has been exercised against real
 |  models and agrees with the arithmetic.  STILL UNEXERCISED:
 |  ETA2_METHOD = DIRECT, MODEL = MIXED, RESIDTYPE = CONDITIONAL and the
 |  degrees-of-freedom diagnostics.  Confirm the output before any of
 |  those goes into a manuscript.
 *===========================================================================*/

%macro mixed_effectsize(tests3    = ,
                        model     = MARGINAL,
                        residtype = ,
                        resids    = ,
                        outpm     = ,
                        dv        = ,
                        resid     = Resid,
                        out       = effectsizes,
                        alpha     = 0.05,
                        label     = ,
                        print     = Y,
                        eta2_method = FROM_F,
                        data      = ,
                        class     = ,
                        fixed     = ,
                        modelinfo = ,
                        dimensions = ,
                        nsubjects = ,
                        between   = );

    %local _abort _havess _haveprobf _cilevel _mdl _rt _src
           _e2 _e2lab _rtword _srcword _e2m _glmok
           _covstruct _ddfm _resvar _subjeff _colz _nsub _havesub _bi;
    %let _abort     = 0;
    %let _glmok     = 0;
    %let _havesub   = 0;
    %let _e2m       = %upcase(%superq(eta2_method));
    %let _covstruct = ;  %let _ddfm = ;  %let _resvar = ;  %let _subjeff = ;
    %let _colz      = ;  %let _nsub = ;
    %let _havess    = 0;
    %let _haveprobf = 0;
    %let _cilevel   = %sysevalf((1 - &alpha) * 100);
    %let _mdl       = %upcase(%superq(model));

    /*-- 0. Resolve the deprecated OUTPM= alias --------------------------*/
    %if %superq(resids) = and %superq(outpm) ne %then %do;
        %let resids = &outpm;
        %put NOTE: (mixed_effectsize) OUTPM= is a deprecated alias for RESIDS= -- using RESIDS=&resids..;
    %end;

    /*-- 1. Validate MODEL= and RESIDTYPE= -------------------------------*/
    %if &_mdl ne MARGINAL and &_mdl ne MIXED %then %do;
        %put ERROR: (mixed_effectsize) MODEL= must be MARGINAL or MIXED, not "&model".;
        %let _abort = 1;
    %end;

    %if %superq(residtype) = %then %let _rt = MARGINAL;
    %else %let _rt = %upcase(%superq(residtype));

    %if &_rt ne MARGINAL and &_rt ne CONDITIONAL %then %do;
        %put ERROR: (mixed_effectsize) RESIDTYPE= must be MARGINAL or CONDITIONAL, not "&residtype".;
        %let _abort = 1;
    %end;

    %if &_mdl = MARGINAL and &_rt = CONDITIONAL %then %do;
        %put ERROR: (mixed_effectsize) RESIDTYPE=CONDITIONAL is meaningless when MODEL=MARGINAL.;
        %put ERROR- A model with no RANDOM statement has no BLUPs, so OUTP= and OUTPM= are identical.;
        %put ERROR- Either drop RESIDTYPE=, or set MODEL=MIXED if the model does have random effects.;
        %let _abort = 1;
    %end;

    /* wording used in notes, titles and labels */
    %if &_rt = MARGINAL %then %do;
        %let _src     = OUTPM=;
        %let _rtword  = marginal;
        %let _srcword = population-average;
        %let _e2      = eta2_approx;
        %let _e2lab   = %str(Eta-squared, proportion of total variance (approximate));
    %end;
    %else %do;
        %let _src     = OUTP=;
        %let _rtword  = conditional;
        %let _srcword = subject-specific;
        %let _e2      = eta2_cond_approx;
        %let _e2lab   = %str(Conditional eta-squared, proportion of WITHIN-SUBJECT variance (approximate));
    %end;

    /*-- 1b. Validate ETA2_METHOD= --------------------------------------*/
    %if &_e2m ne FROM_F and &_e2m ne DIRECT %then %do;
        %put ERROR: (mixed_effectsize) ETA2_METHOD= must be FROM_F or DIRECT, not "&eta2_method".;
        %let _abort = 1;
    %end;

    %if &_e2m = DIRECT %then %do;
        %if %superq(data) = or %superq(fixed) = %then %do;
            %put ERROR: (mixed_effectsize) ETA2_METHOD=DIRECT requires DATA= and FIXED=.;
            %put ERROR- Supply the dataset the model was fitted to and the fixed-effects;
            %put ERROR- right-hand side, exactly as written in the PROC MIXED MODEL statement.;
            %let _abort = 1;
        %end;
        %else %if not %sysfunc(exist(&data)) %then %do;
            %put ERROR: (mixed_effectsize) Dataset &data does not exist.;
            %let _abort = 1;
        %end;
        %if %superq(dv) = %then %do;
            %put ERROR: (mixed_effectsize) ETA2_METHOD=DIRECT requires DV=.;
            %let _abort = 1;
        %end;
        %if %superq(resids) = %then %do;
            %put ERROR: (mixed_effectsize) ETA2_METHOD=DIRECT also requires RESIDS=.;
            %put ERROR- The residuals dataset supplies N_OBS, which the macro compares against the;
            %put ERROR- observation count PROC GLM used and against DenDF. Pass the OUTPM= dataset.;
            %let _abort = 1;
        %end;
        %if &_rt = CONDITIONAL %then %do;
            %put ERROR: (mixed_effectsize) ETA2_METHOD=DIRECT is not defined for RESIDTYPE=CONDITIONAL.;
            %put ERROR- PROC GLM has no random effects, so there is no within-subject decomposition to take.;
            %let _abort = 1;
        %end;
    %end;

    /*-- 2. Validate datasets -------------------------------------------*/
    %if %superq(tests3) = %then %do;
        %put ERROR: (mixed_effectsize) TESTS3= is required.;
        %let _abort = 1;
    %end;
    %else %if not %sysfunc(exist(&tests3)) %then %do;
        %put ERROR: (mixed_effectsize) Dataset &tests3 does not exist.;
        %let _abort = 1;
    %end;

    %if &_abort = 0 %then %do;
        %local dsid rc;
        %let dsid = %sysfunc(open(&tests3));
        %if %sysfunc(varnum(&dsid, NumDF)) = 0 or
            %sysfunc(varnum(&dsid, DenDF)) = 0 or
            %sysfunc(varnum(&dsid, FValue)) = 0 %then %do;
            %put ERROR: (mixed_effectsize) &tests3 must contain NumDF, DenDF and FValue.;
            %put ERROR- Did the PROC MIXED step include  ods output tests3 = &tests3  ?;
            %let _abort = 1;
        %end;
        /* ProbF is optional; note whether it is there so PROC PRINT can adapt */
        %if %sysfunc(varnum(&dsid, ProbF)) > 0 %then %let _haveprobf = 1;
        %let rc = %sysfunc(close(&dsid));
    %end;

    %if %superq(resids) ne %then %do;
        %if %superq(dv) = %then %do;
            %put ERROR: (mixed_effectsize) DV= is required when RESIDS= is supplied.;
            %let _abort = 1;
        %end;
        %else %if not %sysfunc(exist(&resids)) %then %do;
            %put ERROR: (mixed_effectsize) Dataset &resids does not exist.;
            %let _abort = 1;
        %end;
        %else %let _havess = 1;
    %end;

    %if &_abort = 1 %then %do;
        %put ERROR: (mixed_effectsize) Aborting -- no output produced.;
        %return;
    %end;

    /*-- 3. Announce what is being assumed -------------------------------*
     |  The macro cannot tell OUTPM= from OUTP= by inspection, so it states
     |  its assumption loudly rather than guessing.
     *-------------------------------------------------------------------*/
    %put NOTE: (mixed_effectsize) MODEL=&_mdl, RESIDTYPE=&_rt..;
    %if &_havess = 1 %then %do;
        %put NOTE- Assuming &resids contains &_rtword residuals, i.e. the &_src dataset (&_srcword fit).;
        %if &_mdl = MIXED and &_rt = MARGINAL %then
            %put NOTE- For MODEL=MIXED this must be the OUTPM= dataset, NOT OUTP=.;
        %if &_rt = CONDITIONAL %then
            %put NOTE- &_e2 is a WITHIN-SUBJECT proportion and is not comparable with eta-squared from other studies.;
    %end;

    /*-- 3b. Read the model metadata, if it was captured ------------------*
     |  ODS OUTPUT ModelInfo= and Dimensions= cost one line each in the
     |  PROC MIXED step and turn the df checks below from guesswork into
     |  arithmetic.  Everything here is optional; the macro degrades to
     |  its previous behaviour without them.
     *-------------------------------------------------------------------*/
    %if %superq(modelinfo) ne %then %do;
        %if %sysfunc(exist(&modelinfo)) %then %do;
            data _null_;
                set &modelinfo;
                length d $60 v $80;
                d = upcase(strip(vvalue(Descr)));
                v = strip(vvalue(Value));
                if d = 'COVARIANCE STRUCTURE'      then call symputx('_covstruct', v);
                else if d = 'DEGREES OF FREEDOM METHOD' then call symputx('_ddfm', v);
                else if d = 'RESIDUAL VARIANCE METHOD'  then call symputx('_resvar', v);
                else if d = 'SUBJECT EFFECT'           then call symputx('_subjeff', v);
            run;
        %end;
        %else %put WARNING: (mixed_effectsize) MODELINFO=&modelinfo does not exist; df diagnostics limited.;
    %end;

    %if %superq(dimensions) ne %then %do;
        %if %sysfunc(exist(&dimensions)) %then %do;
            data _null_;
                set &dimensions;
                length d $60;
                d = upcase(strip(vvalue(Descr)));
                if d = 'SUBJECTS'      then call symputx('_nsub', input(strip(vvalue(Value)), ?? best32.));
                else if d = 'COLUMNS IN Z' then call symputx('_colz', input(strip(vvalue(Value)), ?? best32.));
            run;
        %end;
        %else %put WARNING: (mixed_effectsize) DIMENSIONS=&dimensions does not exist; df diagnostics limited.;
    %end;

    /* Dimensions "Subjects" is only meaningful when SAS could actually block
       the problem by subject.  A RANDOM effect with no SUBJECT= option spans
       subjects, so SAS reports Subjects=1 and Max Obs per Subject = N -- the
       count is then useless and must not be used as a bound.               */
    %if %superq(nsubjects) ne %then %let _nsub = &nsubjects;   /* explicit always wins */
    %else %if %superq(_nsub) ne %then %do;
        %if %sysevalf(&_nsub <= 1) %then %do;
            %put NOTE: (mixed_effectsize) Dimensions reports Subjects=&_nsub, so SAS could not block;
            %put NOTE- the problem by subject -- usually a RANDOM effect with no SUBJECT= option.;
            %put NOTE- The subject count is not usable as a bound. Pass NSUBJECTS= to enable the;
            %put NOTE- between-subject degrees-of-freedom check.;
            %let _nsub = ;
        %end;
    %end;
    %if %superq(_nsub) ne %then %if %sysevalf(&_nsub > 1) %then %let _havesub = 1;

    %if %superq(_covstruct) ne or &_havesub = 1 %then %do;
        %put NOTE: (mixed_effectsize) Model metadata read:;
        %if %superq(_covstruct) ne %then %put NOTE-   covariance structure  : &_covstruct;
        %if %superq(_resvar)    ne %then %put NOTE-   residual variance     : &_resvar;
        %if %superq(_ddfm)      ne %then %put NOTE-   df method             : &_ddfm;
        %if %superq(_subjeff)   ne %then %put NOTE-   subject effect        : &_subjeff;
        %if %superq(_colz)      ne %then %put NOTE-   columns in Z          : &_colz;
        %if &_havesub = 1          %then %put NOTE-   subjects              : &_nsub;
    %end;

    /*-- 4. Sums of squares from the residuals dataset -------------------*
     |  CSS = corrected SS   = sum((y - ybar)**2)   -> SS_total
     |  USS = uncorrected SS = sum(resid**2)        -> SS_error
     |  USS rather than CSS for the residuals is deliberate: residuals are
     |  centred on zero by construction, and USS is the quantity the
     |  sum-of-squares decomposition calls for.
     *-------------------------------------------------------------------*/
    %if &_havess = 1 %then %do;

        /* WHERE clause is load-bearing.  An OUTPM= dataset holds one row per
           observation the procedure READ, not per observation it USED: rows
           dropped for a missing covariate still appear, with a missing
           residual.  Without this filter CSS(dv) would be taken over a
           different, larger set of rows than USS(resid), and N would count
           rows the model never fitted.  That mismatch -- observations read
           standing in for observations used -- is exactly the error in the
           script this macro replaces.                                      */
        proc means data = &resids noprint nway;
            where not missing(&resid);
            var &dv &resid;
            output out = _es_ss (drop = _type_ _freq_)
                   css(&dv)    = _ss_total
                   uss(&resid) = _ss_error
                   n(&resid)   = _n_obs;
        run;

        /* how many rows the residuals dataset actually holds, for the note below */
        %local _n_rows;
        proc sql noprint;
            select count(*) into :_n_rows trimmed from &resids;
        quit;

        data _null_;
            set _es_ss;
            call symputx('_ss_total', _ss_total);
            call symputx('_ss_error', _ss_error);
            call symputx('_n_obs',    _n_obs);
        run;

        %if %superq(_n_rows) ne and %sysevalf(&_n_rows > &_n_obs) %then %do;
            %put NOTE: (mixed_effectsize) &resids holds &_n_rows rows, of which &_n_obs have a residual.;
            %put NOTE- Sums of squares and N_OBS are taken over the &_n_obs rows the model actually;
            %put NOTE- used. The other rows were read but not fitted -- missing response or;
            %put NOTE- covariate, or timepoint slots belonging to a different outcome in a long file.;
        %end;

        %if %sysevalf(&_ss_total <= 0) %then %do;
            %put WARNING: (mixed_effectsize) SS_total is not positive -- eta-squared suppressed.;
            %let _havess = 0;
        %end;
    %end;

    /*-- 4b. ETA2_METHOD=DIRECT: real Type III SS from PROC GLM ----------*
     |  PROC GLM ignores the repeated-measures covariance structure.  That
     |  is fine here and only here: we take its SUMS OF SQUARES, which are
     |  a genuine decomposition of the observed variance, and we discard
     |  its F and p, which are not valid for this design.  F, df and
     |  partial eta-squared always come from PROC MIXED.
     *-------------------------------------------------------------------*/
    %local _glm_sstot _glm_mse _glm_n;
    %if &_e2m = DIRECT %then %do;

        %put NOTE: (mixed_effectsize) ETA2_METHOD=DIRECT -- refitting the mean model with PROC GLM for Type III SS.;
        %put NOTE- Only the sums of squares are used. GLM F and p values are NOT valid for a repeated-measures design.;

        ods graphics off;
        ods listing close;
        proc glm data = &data;
            %if %superq(class) ne %then %str(class &class;);
            model &dv = &fixed / ss3;
            ods output ModelANOVA = _es_glm  OverallANOVA = _es_glmoa;
        quit;
        ods listing;

        %if %sysfunc(exist(_es_glm)) and %sysfunc(exist(_es_glmoa)) %then %do;

            data _es_glmss;
                set _es_glm;
                where HypothesisType = 3;
                length _eskey $80;
                _eskey  = upcase(compress(Source));
                _ss_eff = SS;
                _df_eff = DF;
                keep _eskey _ss_eff _df_eff;
            run;

            proc sort data = _es_glmss; by _eskey; run;

            data _null_;
                set _es_glmoa;
                if upcase(Source) = 'CORRECTED TOTAL' then do;
                    call symputx('_glm_sstot', SS);
                    call symputx('_glm_n', DF + 1);
                end;
                if upcase(Source) = 'ERROR' then call symputx('_glm_mse', MS);
            run;

            %if %superq(_glm_sstot) ne and %sysevalf(&_glm_sstot > 0) %then %let _glmok = 1;
            %else %put WARNING: (mixed_effectsize) PROC GLM produced no usable corrected total SS; falling back to FROM_F.;

            %if &_glmok = 1 and &_havess = 1 %then %do;
                %if %sysevalf(&_glm_n ne &_n_obs) %then %do;
                    %put WARNING: (mixed_effectsize) PROC GLM used &_glm_n observations, PROC MIXED used &_n_obs..;
                    %put WARNING- The two fits do not have the same complete-case set -- eta-squared and the F test;
                    %put WARNING- then describe slightly different data. Check the missing-value pattern.;
                %end;
            %end;
        %end;
        %else %put WARNING: (mixed_effectsize) PROC GLM did not produce ODS output; falling back to FROM_F.;
    %end;

    %local _e2src;
    %if &_glmok = 1 %then %let _e2src = GLM_TYPE3;
    %else %let _e2src = FROM_F;

    /*-- 5. Effect sizes -------------------------------------------------*/
    data &out;
        length model_label $200 model_type $8 resid_type $11 eta2_source $10 _eskey $80;
        set &tests3;
        model_label = "&label";
        model_type  = "&_mdl";
        resid_type  = "&_rt";
        eta2_source = "&_e2src";
        _eskey      = upcase(compress(Effect));
        _seq        = _n_;   /* preserve Type 3 order across the merge sort */

        %if &_havesub = 1 %then %do;
            /* A BETWEEN-subject effect cannot have more independent units than
               there are subjects, whatever the covariance structure.  This is
               a hard bound; DenDF above it means the test is treating repeated
               observations on the same person as independent.               */
            n_subjects  = &_nsub;
            dendf_ratio = DenDF / &_nsub;
        %end;
        %if %superq(between) ne %then %do;
            length is_between $3;
            is_between = 'no';
            %let _bi = 1;
            %do %while (%scan(%superq(between), &_bi, %str( )) ne );
                if upcase(compress(Effect)) = upcase(compress("%scan(%superq(between), &_bi, %str( ))"))
                    then is_between = 'YES';
                %let _bi = %eval(&_bi + 1);
            %end;
        %end;

        /* ---- partial eta-squared: needs only NumDF, DenDF, FValue ---- */
        if NumDF > 0 and DenDF > 0 and not missing(FValue) then do;
            partial_eta2 = (NumDF * FValue) / (NumDF * FValue + DenDF);
            cohen_f2     = partial_eta2 / (1 - partial_eta2);

            /* Partial epsilon-squared. Like partial eta-squared it
               needs ONLY NumDF, DenDF and FValue, so a reader can
               recompute it from the same table, and it is less
               positively biased. Definition matches R's effectsize
               package (F_to_epsilon2), verified numerically. It can go
               negative when F < 1; that is expected, and such a value
               is conventionally read as 0.                         */
            partial_epsilon2 = (NumDF * (FValue - 1)) / (NumDF * FValue + DenDF);

            /* Steiger (2004) CI by inverting the noncentral F.
               FNONCT(x, ndf, ddf, p) returns the noncentrality nc for which
               PROBF(x, ndf, ddf, nc) = p.                                */
            _n_eff = NumDF + DenDF + 1;

            if probf(FValue, NumDF, DenDF, 0) > (1 - &alpha/2)
                then _nc_lo = fnonct(FValue, NumDF, DenDF, 1 - &alpha/2);
                else _nc_lo = 0;
            _nc_hi = fnonct(FValue, NumDF, DenDF, &alpha/2);

            partial_eta2_lcl = _nc_lo / (_nc_lo + _n_eff);
            partial_eta2_ucl = _nc_hi / (_nc_hi + _n_eff);
        end;

        /* ---- eta-squared: needs the residuals ------------------------ */
        %if &_havess = 1 %then %do;
            ss_total = &_ss_total;
            ss_error = &_ss_error;
            n_obs    = &_n_obs;

            if NumDF > 0 and DenDF > 0 and not missing(FValue) then do;
                mse       = ss_error / DenDF;
                ss_effect = NumDF * FValue * mse;

                &_e2 = ss_effect / ss_total;

                /* How far is DenDF from the observation count?  See
                   CHOICE 3 in the header: it makes the MSE substitution
                   visible but does NOT bound the eta-squared error.
                   Partial eta-squared is immune either way.          */
                nobs_dendf = n_obs / DenDF;
            end;
        %end;

        drop _nc_lo _nc_hi _n_eff;

        label
            Effect           = 'Effect'
            NumDF            = 'Numerator df'
            DenDF            = 'Denominator df'
            FValue           = 'F value'
            ProbF            = 'Pr > F'
            partial_eta2     = 'Partial eta-squared'
            partial_eta2_lcl = 'Partial eta-squared, lower CL'
            partial_eta2_ucl = 'Partial eta-squared, upper CL'
            cohen_f2         = "Cohen's f-squared"
            partial_epsilon2 = 'Partial epsilon-squared (less biased; from F and df only)'
        %if &_havess = 1 %then %do;
            &_e2             = "&_e2lab"
            mse              = 'Implied MSE (SS_error / DenDF)'
            ss_effect        = 'SS effect (approximate)'
            ss_error         = "SS error (USS of &_rtword residuals)"
            ss_total         = 'SS total (CSS of dependent variable)'
            n_obs            = 'Observations used'
            nobs_dendf       = 'N_obs / DenDF (>1.5 : eta-squared from F unreliable)'
        %end;
        %if &_havesub = 1 %then %do;
            n_subjects       = 'Subjects'
            dendf_ratio      = 'DenDF / subjects (>1 : more df than independent units)'
        %end;
        %if %superq(between) ne %then %do;
            is_between       = 'Between-subject effect?'
        %end;
            eta2_source      = 'Source of eta-squared'
            model_label      = 'Model'
            model_type       = 'Model type'
            resid_type       = 'Residual type'
        ;
    run;

    /*-- 5b. Overwrite eta-squared with the DIRECT Type III values -------*/
    %if &_glmok = 1 %then %do;

        proc sort data = &out; by _eskey; run;

        data &out;
            merge &out (in = a) _es_glmss (in = b);
            by _eskey;
            if a;
            if b then do;
                ss_effect  = _ss_eff;
                ss_total   = &_glm_sstot;
                mse        = &_glm_mse;
                &_e2 = ss_effect / ss_total;
            end;
            else do;
                call missing(ss_effect, &_e2);
                eta2_source = 'UNMATCHED';
                put "WARNING: (mixed_effectsize) No Type III SS matched effect " Effect
                    "-- eta-squared suppressed for that row.";
            end;
            label mse       = 'MSE (PROC GLM Type III error)'
                  ss_effect = 'SS effect (PROC GLM Type III)'
                  ss_total  = 'SS total (PROC GLM corrected total)';
            drop _ss_eff _df_eff;
        run;

        %put NOTE: (mixed_effectsize) Eta-squared taken from PROC GLM Type III sums of squares.;
        %put NOTE- These are a real variance decomposition, not reconstructed from F. Partial eta-squared is unchanged.;
    %end;

    /* restore the original Type 3 ordering and drop the merge keys */
    proc sort data = &out; by _seq; run;
    data &out;
        set &out;
        drop _eskey _seq;
    run;

    /*-- 6a. Warn when eta-squared from F cannot be trusted ---------------*
     |  DenDF much smaller than N_obs means the implied MSE (SS_error /
     |  DenDF) is inflated, and eta-squared with it.  That is the normal
     |  situation for a BETWEEN-subject effect in a repeated-measures
     |  design.  Partial eta-squared is unaffected.
     *-------------------------------------------------------------------*/
    %local _worst;
    %if &_havess = 1 and &_glmok = 0 %then %do;
        %let _worst = .;
        proc sql noprint;
            select max(nobs_dendf) into :_worst trimmed
              from &out where nobs_dendf is not missing;
        quit;
        %put WARNING: (mixed_effectsize) ETA2_METHOD=FROM_F reconstructs SS_effect as NumDF x F x MSE,;
        %put WARNING- which is exact only in a fixed-effects ANOVA. Here it can be out by several fold;
        %put WARNING- in EITHER direction -- worst N_obs/DenDF on this table: &_worst -- and the error;
        %put WARNING- is NOT predictable from the df. Report PARTIAL_ETA2 or PARTIAL_EPSILON2,;
        %put WARNING- which are unaffected, or rerun with ETA2_METHOD=DIRECT.;
        %put WARNING- The measured errors: README.md, "Caveats to state in a Methods section", item 2.;
    %end;

    /*-- 6ab. Between-subject effects tested on observation-scale df -------*
     |  A BETWEEN-subject effect has at most one independent unit per
     |  subject, so a DenDF above the subject count puts it on the
     |  OBSERVATION scale.  The TEST is not wrong -- the model-based SE
     |  already uses the fitted covariance -- but PARTIAL_ETA2 is then
     |  conditioned on the df method.  See README.md, "Choose a
     |  denominator-df method", for what to report.
     *-------------------------------------------------------------------*/
    %if &_havesub = 1 %then %do;
        %local _nover _worstratio _overlist;
        %let _nover = 0; %let _worstratio = .; %let _overlist = ;
        proc sql noprint;
            select count(*), max(dendf_ratio) into :_nover trimmed, :_worstratio trimmed
              from &out where dendf_ratio > 1
              %if %superq(between) ne %then %str(and is_between = 'YES');
              ;
            select Effect into :_overlist separated by ', '
              from &out where dendf_ratio > 1
              %if %superq(between) ne %then %str(and is_between = 'YES');
              ;
        quit;
        %if %superq(_nover) ne and %superq(_nover) ne . %then %do;
        %if &_nover > 0 %then %do;
            %put WARNING: (mixed_effectsize) &_nover effect(s) have DenDF above the subject count (&_nsub).;
            %put WARNING- Effects: &_overlist;
            %put WARNING- Worst DenDF/subjects ratio: &_worstratio..;
            %if %superq(between) ne %then %do;
            %put WARNING- Declared BETWEEN-subject via BETWEEN=, so these are on OBSERVATION-scale df.;
            %end;
            %else %do;
            %put WARNING- Pass BETWEEN= to sharpen this. A WITHIN-subject effect may legitimately;
            %put WARNING- exceed the subject count under compound symmetry, but not under UN.;
            %end;
            %put WARNING- The TEST is unaffected -- the GLS standard error already uses the fitted;
            %put WARNING- covariance. PARTIAL_ETA2 is CONDITIONED on the df method and has no df-free;
            %put WARNING- value: state the method beside it, or report an estimate with a CI instead.;
            %put WARNING- Why: README.md, "Choose a denominator-df method".;
            %if %superq(_ddfm) ne %then
            %put WARNING- Degrees of Freedom Method in this fit: &_ddfm..;
        %end;
        %end;
    %end;

    /*-- 6ac. Containment fall-through -------------------------------------*/
    %if %superq(_ddfm) ne %then %do;
        %if %index(%upcase(&_ddfm), CONTAIN) > 0 %then %do;
            %put NOTE: (mixed_effectsize) DDFM=CONTAINMENT is in force -- the SAS default whenever a;
            %put NOTE- RANDOM statement is present. An effect contained by no RANDOM effect falls;
            %put NOTE- through to the residual df, N - rank(XZ), on the OBSERVATION scale, which a;
            %put NOTE- subject-level predictor usually does. LaBGAS standard is DDFM=KR2 -- README.md.;
            %if %superq(_colz) ne %then %put NOTE-   Columns in Z = &_colz (a RANDOM statement is active).;
        %end;
    %end;

    /*-- 6b. Warn if the denominator df look like a default ---------------*/
    %local _fracdf;
    %let _fracdf = .;
    proc sql noprint;
        select sum(DenDF ne int(DenDF)) into :_fracdf trimmed from &out;
    quit;
    %if %superq(_fracdf) ne and %superq(_fracdf) ne . %then %do;
        %if &_fracdf = 0 %then %do;
            %put NOTE: (mixed_effectsize) Every DenDF is a whole number, so the df were probably left;
            %put NOTE- at the SAS default -- neither KR/KR2 nor SATTERTH was requested. The LaBGAS;
            %put NOTE- standard is DDFM=KR2 on every mixed or marginal model. If it was deliberately;
            %put NOTE- not used, record why. Under type=UN the default df are close to exact and need;
            %put NOTE- no correction -- see README.md, "Choose a denominator-df method".;
        %end;
    %end;

    /* One DenDF on every row: EXPECTED under UN, suspicious under CS/VC/AR(1).
       The Type 3 table does not carry the covariance structure, so the macro
       cannot tell which case it is in and must not call this an error.      */
    %local _ndistinct _neff;
    proc sql noprint;
        select count(distinct DenDF), count(*) into :_ndistinct trimmed, :_neff trimmed from &out;
    quit;
    %if %superq(_ndistinct) ne and %superq(_neff) ne %then %do;
        %if &_ndistinct = 1 and &_neff > 1 %then %do;
            %put NOTE: (mixed_effectsize) All &_neff effects share the same DenDF.;
            %put NOTE- Under an UNSTRUCTURED covariance that is EXPECTED and very nearly correct --;
            %put NOTE- do not "fix" it. Under CS, VC or AR(1) WITH a residual term it does indicate;
            %put NOTE- that the df method is not stratifying. Check Model Information, and see;
            %put NOTE- README.md, "Choose a denominator-df method".;
        %end;
    %end;

    /*-- 7. Print --------------------------------------------------------*/
    %if %upcase(&print) = Y %then %do;
        title  "Effect sizes for fixed effects%if %superq(label) ne %then %str( -- &label);";
        title2 "&_mdl model; partial eta-squared = (NumDF x F) / (NumDF x F + DenDF), with &_cilevel%str(%%) CI by noncentral F";
        %if &_havess = 1 %then %do;
        title3 "&_e2 computed from &_rtword residuals -- an approximation for a mixed model; see README.md";
        %end;

        proc print data = &out noobs label;
            var Effect NumDF DenDF FValue
                %if &_haveprobf = 1 %then ProbF;
                partial_eta2 partial_eta2_lcl partial_eta2_ucl
                partial_epsilon2 cohen_f2
                %if &_havess = 1 %then &_e2 nobs_dendf;
                %if &_havesub = 1 %then dendf_ratio;
                %if %superq(between) ne %then is_between;
                eta2_source
                ;
            format FValue 8.2 NumDF 8. DenDF 8.1
                   %if &_haveprobf = 1 %then %str(ProbF pvalue6.4);
                   partial_eta2 partial_eta2_lcl partial_eta2_ucl
                   partial_epsilon2 cohen_f2 8.3
                   %if &_havess = 1 %then %str(&_e2 8.4 nobs_dendf 8.2);
                   %if &_havesub = 1 %then %str(dendf_ratio 8.2);
                   ;
        run;
        title;
    %end;

    proc datasets library = work nolist nowarn;
        delete _es_ss _es_glm _es_glmoa _es_glmss;
    quit;

    %put NOTE: (mixed_effectsize) Effect sizes written to &out..;
    %put NOTE- Report PARTIAL_ETA2, optionally with PARTIAL_EPSILON2, which is less biased and is;
    %put NOTE- recomputable by a reader from the F and df in the same table.;
    %put NOTE- If you also report &_e2, label it as the output labels do -- it is the CLASSICAL,;
    %put NOTE- non-partial quantity and is not comparable with PARTIAL_ETA2.;

%mend mixed_effectsize;
