/*===========================================================================*
 |  mixed_effectsize.sas
 |
 |  Effect sizes for fixed effects in models fitted with PROC MIXED --
 |  both MARGINAL models (REPEATED statement, population-average) and true
 |  MIXED models (RANDOM statement, subject-specific).
 |
 |  Produces, for every effect in the Type 3 tests table:
 |      partial_eta2   partial eta-squared, with a confidence interval
 |      cohen_f2       Cohen's f-squared (a monotone function of the above)
 |      eta2_*         eta-squared     (optional; needs a residuals dataset)
 |      omega2_*       omega-squared   (optional; needs a residuals dataset)
 |
 |  Every quantity is written to the output dataset under its own name and
 |  printed with a full label, so that nothing can be reported under the
 |  wrong heading.
 |
 |  Companion file:  es_identify.sas  -- forensic helper for auditing effect
 |  sizes in tables that were produced by other code.
 |
 |  ---------------------------------------------------------------------
 |  WHY THIS EXISTS
 |  ---------------------------------------------------------------------
 |  An earlier in-house script computed three statistics and it was easy to
 |  report one of them under the name of another. It also contained a bug:
 |  the error sum of squares was formed as
 |
 |      ss_error = mse * (_freq_ - numdf);
 |
 |  which uses the NUMBER OF OBSERVATIONS minus the effect's NUMERATOR
 |  degrees of freedom, where the model's DENOMINATOR degrees of freedom
 |  belong.  In a repeated-measures design with ~660 observations and a
 |  denominator df of ~161, that inflates the denominator roughly fourfold
 |  and shrinks partial eta-squared by the same factor.
 |
 |  This macro never touches _FREQ_.  Partial eta-squared needs nothing but
 |  the Type 3 table, because PROC MIXED already reports the Kenward-Roger
 |  or Satterthwaite denominator df there:
 |
 |      partial_eta2 = (NumDF * FValue) / (NumDF * FValue + DenDF)
 |
 |  This is the standard partial eta-squared identity.  It is also the
 |  quantity Edwards, Muller, Wolfinger, Qaqish & Schabenberger (2008)
 |  propose as R-squared-beta for fixed effects in the linear mixed model
 |  (Stat Med 27:6137-6157, doi:10.1002/sim.3429).  Check the exact
 |  correspondence in that paper before citing it in a Methods section.
 |
 |  ---------------------------------------------------------------------
 |  MARGINAL versus MIXED -- what actually changes
 |  ---------------------------------------------------------------------
 |  PARTIAL ETA-SQUARED IS THE SAME IN BOTH CASES.  The formula depends
 |  only on NumDF, DenDF and F, and the denominator df already reflects
 |  whatever covariance structure was fitted.  Nothing to choose.
 |
 |  ETA-SQUARED AND OMEGA-SQUARED DO CHANGE, because "the residual" means
 |  two different things:
 |
 |    MARGINAL residuals   y - Xb          from the OUTPM= dataset.
 |                         Variation around the population-average fit.
 |                         SS_error still CONTAINS between-subject
 |                         variance.  Eta-squared is then the proportion
 |                         of TOTAL observed variance, which is the honest
 |                         analogue of classical eta-squared and is
 |                         comparable with fixed-effects studies.
 |
 |    CONDITIONAL residuals  y - Xb - Zu   from the OUTP= dataset.
 |                         Variation around the subject-specific fit.
 |                         The random effects have already absorbed the
 |                         between-subject variance, so SS_error is much
 |                         smaller and the resulting quantity is a
 |                         WITHIN-SUBJECT proportion of variance.  It is
 |                         systematically larger and is NOT comparable
 |                         with eta-squared from other studies.
 |
 |  The macro names the outputs differently in the two cases
 |  (ETA2_APPROX versus ETA2_COND_APPROX) precisely so that the two cannot
 |  be confused in a results table.
 |
 |  DEFAULT: marginal residuals, for both model types.  That is the
 |  comparable quantity.  Ask for conditional residuals only when the
 |  within-subject proportion is what you actually want to report, and say
 |  so explicitly in the Methods.
 |
 |  For a marginal model there are no random effects, so OUTP= and OUTPM=
 |  are identical and RESIDTYPE=CONDITIONAL is refused as meaningless.
 |
 |  ---------------------------------------------------------------------
 |  USAGE
 |  ---------------------------------------------------------------------
 |  Marginal model (REPEATED):
 |
 |      proc mixed data = mydata;
 |          class id group time;
 |          model y = group|time / outpm = outpm_y ddfm = kr;
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
 |          model y = group|time / outpm = outpm_y outp = outp_y ddfm = kr;
 |          random intercept / subject = id;
 |          ods output tests3 = tests3_y;
 |      run;
 |
 |      %mixed_effectsize(tests3 = tests3_y,
 |                        model  = MIXED,
 |                        resids = outpm_y,        /* note: OUTPM, not OUTP */
 |                        dv     = y);
 |
 |  True mixed model, within-subject (conditional) eta-squared:
 |
 |      %mixed_effectsize(tests3    = tests3_y,
 |                        model     = MIXED,
 |                        residtype = CONDITIONAL,
 |                        resids    = outp_y,      /* note: OUTP */
 |                        dv        = y);
 |
 |  Partial eta-squared and its CI need only TESTS3.  RESIDS= and DV= are
 |  needed only for eta-squared and omega-squared.
 |
 |  ---------------------------------------------------------------------
 |  PARAMETERS
 |  ---------------------------------------------------------------------
 |  tests3     Dataset from ODS OUTPUT TESTS3=.  Required.
 |             Must contain Effect, NumDF, DenDF, FValue (ProbF optional).
 |  model      MARGINAL (default) for a REPEATED-only model, MIXED for a
 |             model with a RANDOM statement.  Controls validation and the
 |             wording of the output; the partial eta-squared arithmetic is
 |             identical either way.
 |  residtype  MARGINAL (default) or CONDITIONAL.  Determines which
 |             residuals the supplied dataset is assumed to contain, and
 |             hence how eta-squared is named and interpreted.
 |             CONDITIONAL is refused when MODEL = MARGINAL.
 |  resids     Residuals dataset: the OUTPM= dataset when
 |             RESIDTYPE = MARGINAL, the OUTP= dataset when CONDITIONAL.
 |             Optional; omit if only partial eta-squared is wanted.
 |  outpm      Deprecated alias for RESIDS=, kept so older calls still run.
 |  dv         Name of the dependent variable in the residuals dataset.
 |             Required if RESIDS= is supplied.
 |  resid      Name of the residual variable.  Default Resid (correct for
 |             both OUTPM= and OUTP=).
 |  out        Output dataset.  Default EFFECTSIZES.
 |  alpha      Two-sided alpha for the partial eta-squared CI.  Default 0.05.
 |  label      Free-text model label carried into the output.  Optional.
 |  print      Y (default) prints a formatted table; N suppresses it.
 |
 |  ---------------------------------------------------------------------
 |  WHICH ONE TO REPORT
 |  ---------------------------------------------------------------------
 |  Report PARTIAL_ETA2.  It is the convention in this literature, and a
 |  reader can recompute it from the F and df printed in the same table,
 |  which is exactly what makes a results table checkable.
 |
 |  Eta-squared cannot be recovered from F and df alone, so a reader cannot
 |  verify it.  Report it only alongside partial eta-squared, never instead
 |  of it, and name it correctly.
 |
 |  Consider also reporting a standardised mean difference with a CI for
 |  the group contrasts specifically -- ESTIMATE or LSMESTIMATE with a
 |  sensible divisor.  For a two-group comparison that is more directly
 |  interpretable than any variance-explained measure.
 |
 |  ---------------------------------------------------------------------
 |  ASSUMPTIONS AND LIMITS  (state these in a Methods section)
 |  ---------------------------------------------------------------------
 |  1. Eta-squared and omega-squared are NOT uniquely defined for a mixed
 |     or marginal model -- there is no single residual variance and no
 |     single error df.  Here they are formed as
 |
 |         SS_error  = USS(residuals)                from RESIDS=
 |         SS_total  = CSS(dependent variable)       from RESIDS=
 |         MSE       = SS_error / DenDF
 |         SS_effect = NumDF * FValue * MSE
 |
 |     which is the natural analogue of the fixed-effects definitions but
 |     is an approximation.  Both are labelled *_APPROX for that reason.
 |
 |  2. Under Kenward-Roger, DenDF differs between effects, so the implied
 |     MSE differs slightly between effects too.  Partial eta-squared is
 |     unaffected (the MSE cancels).  Eta-squared and omega-squared are
 |     mildly affected; the macro reports the DenDF used for each row so
 |     the dependence is visible.
 |
 |  3. With MODEL = MIXED and RESIDTYPE = CONDITIONAL, report the variance
 |     components alongside the effect sizes.  Without them a reader cannot
 |     tell how much variance the random effects absorbed, and the
 |     within-subject proportion is uninterpretable on its own.
 |
 |  4. The CI for partial eta-squared inverts the noncentral F
 |     distribution (Steiger 2004) and uses N = NumDF + DenDF + 1.  With
 |     Kenward-Roger df that N is an approximation, so treat the interval
 |     as approximate.
 |
 |  5. Use DDFM = KR (or KR2) in the MODEL statement.  With the default
 |     containment method the denominator df can be badly wrong for these
 |     designs, and every quantity here depends on it.  The macro prints a
 |     NOTE if every DenDF is a whole number, which usually means KR or
 |     Satterthwaite was not requested.
 |
 |  ---------------------------------------------------------------------
 |  CORRESPONDENCE WITH PROC GLM
 |  ---------------------------------------------------------------------
 |  PROC GLM's EFFECTSIZE option defines the partial correlation ratio as
 |
 |      partial eta-squared = SS_effect / (SS_effect + SS_error)
 |
 |  The formula used here is algebraically identical to it.  In GLM,
 |  F = (SS_effect/df_effect) / MSE with MSE = SS_error/df_error, so
 |  SS_effect = df_effect*F*MSE and SS_error = df_error*MSE.  Substituting,
 |  the MSE cancels:
 |
 |      SS_effect/(SS_effect + SS_error)
 |          = (df_effect*F*MSE) / (df_effect*F*MSE + df_error*MSE)
 |          = (NumDF*F) / (NumDF*F + DenDF)
 |
 |  Verified numerically on an unbalanced 2x3 factorial with Type III sums
 |  of squares: the two expressions agree to machine precision (largest
 |  discrepancy 5.6e-17 across three effects).  The semipartial (plain)
 |  eta-squared, SS_effect/SS_total, likewise agrees exactly, which
 |  confirms that reconstructing SS_effect as NumDF*FValue*MSE is exact in
 |  the fixed-effects case.
 |
 |  TWO THINGS THAT "IDENTICAL FORMULA" DOES NOT MEAN:
 |
 |    1. In GLM, DenDF is the single residual df, the same for every
 |       effect.  Under Kenward-Roger in MIXED it is effect-specific and
 |       usually fractional.  That is correct -- it is the right df for
 |       that effect's test -- but the values are not slices of one sum-of-
 |       squares decomposition the way GLM's are.
 |
 |    2. In a mixed or marginal model there is NO exact sum-of-squares
 |       decomposition at all.  With correlated errors the F test is a
 |       Wald-type test on generalised least squares estimates, not a ratio
 |       of orthogonal sums of squares.  So the identity runs one way only:
 |       the formula can be computed, but the result cannot be recovered as
 |       SS_effect/(SS_effect + SS_error) from any real decomposition.
 |       Read it as the effect size implied by the F test -- which is
 |       precisely what Edwards et al. (2008) formalise as R-squared-beta.
 |       Worth one sentence in a Methods section.
 |
 |  ---------------------------------------------------------------------
 |  OPEN QUESTION: WHICH OMEGA-SQUARED CONVENTION
 |  ---------------------------------------------------------------------
 |  The omega-squared computed here uses
 |
 |      (SS_effect - NumDF*MSE) / (SS_total + MSE)
 |
 |  inherited from the 2021 script.  That is the classic Hays / Keppel
 |  convention, and is also what Olejnik & Algina (2003) use.
 |
 |  The SAS documentation for the EFFECTSIZE option renders semipartial
 |  omega-squared with SS_total alone in the denominator in one place, and
 |  describes it as "the total sum of squares adjusted by MSE terms" in
 |  another.  The formulas on those pages are published as images rather
 |  than text, so THIS HAS NOT BEEN VERIFIED against SAS's actual output.
 |
 |  The difference is small -- about 0.8 per cent in a test case -- and partial
 |  eta-squared is unaffected either way.  To settle it, run PROC GLM with
 |  the EFFECTSIZE option on any small dataset and compare the printed
 |  semipartial omega-squared against both:
 |
 |      (SS_effect - df_effect*MSE) / SS_total
 |      (SS_effect - df_effect*MSE) / (SS_total + MSE)
 |
 |  then update this header and the README with what you find.
 |
 |  ---------------------------------------------------------------------
 |  REFERENCES
 |  ---------------------------------------------------------------------
 |  Edwards LJ, Muller KE, Wolfinger RD, Qaqish BF, Schabenberger O.
 |      An R2 statistic for fixed effects in the linear mixed model.
 |      Statistics in Medicine 2008;27(29):6137-6157.
 |      doi:10.1002/sim.3429   (PMID 18816511)
 |
 |  Olejnik S, Algina J.  Generalized eta and omega squared statistics:
 |      measures of effect size for some common research designs.
 |      Psychological Methods 2003;8(4):434-447.
 |      doi:10.1037/1082-989X.8.4.434
 |
 |  Steiger JH.  Beyond the F test: effect size confidence intervals and
 |      tests of close fit in the analysis of variance and contrast
 |      analysis.  Psychological Methods 2004;9(2):164-182.
 |      doi:10.1037/1082-989X.9.2.164   (PMID 15137887)
 |
 |  SAS Institute.  Effect Size Measures for F Tests in GLM.
 |      SAS/STAT User's Guide.
 |      https://support.sas.com/documentation/cdl/en/statug/63962/HTML/default/statug_glm_sect032.htm
 |
 |  ---------------------------------------------------------------------
 |  Written for the Laboratory for Brain-Gut Axis Studies (LaBGAS),
 |  KU Leuven.  Supersedes "Calculation of effect size following marginal
 |  linear mixed models" (B. Dalile, 3 December 2021).
 |
 |  NOT YET RUN IN SAS -- the arithmetic has been verified independently,
 |  but please run it against a known model and confirm the output before
 |  using it in a manuscript.
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
                        print     = Y);

    %local _abort _havess _haveprobf _cilevel _mdl _rt _src
           _e2 _o2 _e2lab _o2lab _rtword _srcword;
    %let _abort     = 0;
    %let _havess    = 0;
    %let _haveprobf = 0;
    %let _cilevel   = %sysevalf((1 - &alpha) * 100);
    %let _mdl       = %upcase(%superq(model));

    /*-- 0. Resolve the deprecated OUTPM= alias --------------------------*/
    %if %superq(resids) = and %superq(outpm) ne %then %do;
        %let resids = &outpm;
        %put NOTE: (mixed_effectsize) OUTPM= is a deprecated alias for RESIDS=; using RESIDS=&resids..;
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
        %let _o2      = omega2_approx;
        %let _e2lab   = %str(Eta-squared, proportion of total variance (approximate));
        %let _o2lab   = %str(Omega-squared, proportion of total variance (approximate));
    %end;
    %else %do;
        %let _src     = OUTP=;
        %let _rtword  = conditional;
        %let _srcword = subject-specific;
        %let _e2      = eta2_cond_approx;
        %let _o2      = omega2_cond_approx;
        %let _e2lab   = %str(Conditional eta-squared, proportion of WITHIN-SUBJECT variance (approximate));
        %let _o2lab   = %str(Conditional omega-squared, proportion of WITHIN-SUBJECT variance (approximate));
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
            %put ERROR- Did you write ODS OUTPUT TESTS3=&tests3; in the PROC MIXED step?;
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
        %put ERROR: (mixed_effectsize) Aborting; no output produced.;
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

    /*-- 4. Sums of squares from the residuals dataset -------------------*
     |  CSS = corrected SS   = sum((y - ybar)**2)   -> SS_total
     |  USS = uncorrected SS = sum(resid**2)        -> SS_error
     |  USS rather than CSS for the residuals is deliberate: residuals are
     |  centred on zero by construction, and USS is the quantity the
     |  sum-of-squares decomposition calls for.
     *-------------------------------------------------------------------*/
    %if &_havess = 1 %then %do;

        proc means data = &resids noprint nway;
            var &dv &resid;
            output out = _es_ss (drop = _type_ _freq_)
                   css(&dv)    = _ss_total
                   uss(&resid) = _ss_error
                   n(&dv)      = _n_obs;
        run;

        data _null_;
            set _es_ss;
            call symputx('_ss_total', _ss_total);
            call symputx('_ss_error', _ss_error);
            call symputx('_n_obs',    _n_obs);
        run;

        %if %sysevalf(&_ss_total <= 0) %then %do;
            %put WARNING: (mixed_effectsize) SS_total is not positive; eta-squared suppressed.;
            %let _havess = 0;
        %end;
    %end;

    /*-- 5. Effect sizes -------------------------------------------------*/
    data &out;
        length model_label $200 model_type $8 resid_type $11;
        set &tests3;
        model_label = "&label";
        model_type  = "&_mdl";
        resid_type  = "&_rt";

        /* ---- partial eta-squared: needs only NumDF, DenDF, FValue ---- */
        if NumDF > 0 and DenDF > 0 and not missing(FValue) then do;
            partial_eta2 = (NumDF * FValue) / (NumDF * FValue + DenDF);
            cohen_f2     = partial_eta2 / (1 - partial_eta2);

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

        /* ---- eta-squared and omega-squared: need the residuals ------- */
        %if &_havess = 1 %then %do;
            ss_total = &_ss_total;
            ss_error = &_ss_error;
            n_obs    = &_n_obs;

            if NumDF > 0 and DenDF > 0 and not missing(FValue) then do;
                mse       = ss_error / DenDF;
                ss_effect = NumDF * FValue * mse;

                &_e2 = ss_effect / ss_total;
                &_o2 = (ss_effect - NumDF * mse) / (ss_total + mse);
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
        %if &_havess = 1 %then %do;
            &_e2             = "&_e2lab"
            &_o2             = "&_o2lab"
            mse              = 'Implied MSE (SS_error / DenDF)'
            ss_effect        = 'SS effect (approximate)'
            ss_error         = "SS error (USS of &_rtword residuals)"
            ss_total         = 'SS total (CSS of dependent variable)'
            n_obs            = 'Observations used'
        %end;
            model_label      = 'Model'
            model_type       = 'Model type'
            resid_type       = 'Residual type'
        ;
    run;

    /*-- 6. Warn if the denominator df look like the containment default -*/
    %local _fracdf;
    %let _fracdf = .;
    proc sql noprint;
        select sum(DenDF ne int(DenDF)) into :_fracdf trimmed from &out;
    quit;
    %if %superq(_fracdf) ne and %superq(_fracdf) ne . %then %do;
        %if &_fracdf = 0 %then %do;
            %put WARNING: (mixed_effectsize) Every DenDF is a whole number.;
            %put WARNING- That usually means DDFM=KR or DDFM=SATTERTH was not requested.;
            %put WARNING- All effect sizes here depend on DenDF; consider refitting with DDFM=KR.;
        %end;
    %end;

    /*-- 7. Print --------------------------------------------------------*/
    %if %upcase(&print) = Y %then %do;
        title  "Effect sizes for fixed effects%if %superq(label) ne %then %str( -- &label);";
        title2 "&_mdl model; partial eta-squared = (NumDF x F) / (NumDF x F + DenDF), with &_cilevel%str(%%) CI by noncentral F";
        %if &_havess = 1 %then %do;
        title3 "&_e2 computed from &_rtword residuals -- an approximation for a mixed model; see the macro header";
        %end;

        proc print data = &out noobs label;
            var Effect NumDF DenDF FValue
                %if &_haveprobf = 1 %then ProbF;
                partial_eta2 partial_eta2_lcl partial_eta2_ucl cohen_f2
                %if &_havess = 1 %then &_e2 &_o2;
                ;
            format FValue 8.2 NumDF 8. DenDF 8.1
                   %if &_haveprobf = 1 %then %str(ProbF pvalue6.4);
                   partial_eta2 partial_eta2_lcl partial_eta2_ucl
                   cohen_f2 8.3
                   %if &_havess = 1 %then %str(&_e2 &_o2 8.4);
                   ;
        run;
        title;
    %end;

    proc datasets library = work nolist nowarn;
        delete _es_ss;
    quit;

    %put NOTE: (mixed_effectsize) Effect sizes written to &out..;
    %put NOTE- Report PARTIAL_ETA2. If you also report &_e2, label it as the output labels do.;

%mend mixed_effectsize;
