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
 |  --- eta-squared method (does not affect partial eta-squared) --------
 |  eta2_method  FROM_F (default) reconstructs SS_effect from the F test,
 |             as NumDF * FValue * (SS_error / DenDF).  Cheap, but see
 |             ASSUMPTIONS 1-2: it is unreliable when DenDF is much
 |             smaller than the number of observations, which is exactly
 |             the case for a BETWEEN-subject effect in a repeated-
 |             measures design.
 |             DIRECT instead refits the same mean model with PROC GLM
 |             and takes real Type III sums of squares.  Slower, but it
 |             is an actual variance decomposition rather than one
 |             back-solved from a Wald statistic.  Requires DATA=,
 |             CLASS= and FIXED=.  Partial eta-squared is unchanged.
 |  data       Input dataset for ETA2_METHOD = DIRECT.  Use the same
 |             dataset the PROC MIXED model was fitted to.
 |  class      CLASS variables for the DIRECT refit, space separated.
 |  fixed      Right-hand side of the MODEL statement for the DIRECT
 |             refit -- the fixed effects only, exactly as written in
 |             PROC MIXED.
 |
 |  NOTE on DIRECT: PROC GLM ignores the repeated-measures covariance
 |  structure entirely.  Only its SUMS OF SQUARES are used; its F values
 |  and p values are NOT valid for these designs and are never reported
 |  by this macro.  The F, df and partial eta-squared always come from
 |  PROC MIXED.
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
 |  2. *** ETA-SQUARED FROM F IS UNRELIABLE WHEN DenDF << N_obs. ***
 |     MSE is formed as SS_error / DenDF, where SS_error is the USS of
 |     the residuals over ALL observations.  In a fixed-effects ANOVA
 |     that is exact, because the error df there IS N - rank(X).  In a
 |     marginal or mixed model it is not: a BETWEEN-subject effect in a
 |     repeated-measures design is tested against participant-level df,
 |     so DenDF can be several times smaller than N_obs, and the implied
 |     MSE -- and hence eta-squared and omega-squared -- is inflated by
 |     roughly N_obs / DenDF.
 |
 |     Worked case: 166 participants x 4 timepoints = 664 observations,
 |     group effect tested on 161 df.  Eta-squared from F is too large by
 |     a factor of about 664/161 = 4.1.
 |
 |     The macro computes N_OBS / DenDF for every row, reports it as
 |     NOBS_DENDF, and issues a WARNING for any row above 1.5.  When that
 |     warning fires, report PARTIAL_ETA2 only, or rerun with
 |     ETA2_METHOD = DIRECT.
 |
 |     PARTIAL_ETA2 IS NOT AFFECTED -- the MSE cancels out of it
 |     algebraically.  This caveat is about eta-squared and omega-squared
 |     alone.
 |
 |     Under Kenward-Roger, DenDF also differs from effect to effect, so
 |     the implied MSE differs between rows of the same model.  The macro
 |     reports the DenDF used for each row so the dependence is visible.
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
 |  5. Choose a denominator-df method deliberately; do not leave it at
 |     the default for a between-within design.  Everything here depends
 |     on DenDF, so this is the single most consequential choice upstream
 |     of the macro.
 |
 |     The point is NOT that Kenward-Roger is universally best.  It is
 |     that the default containment / between-within methods do not
 |     account for the fact that the covariance parameters were
 |     ESTIMATED, and can be markedly anticonservative in small,
 |     unbalanced designs with a rich covariance structure.
 |
 |       DDFM = KR    Kenward-Roger.  Does two things: inflates the
 |                    covariance matrix of the fixed effects to allow for
 |                    estimated variance parameters, AND derives a
 |                    Satterthwaite-type df from it.  So it changes both
 |                    F and DenDF.  Designed for METHOD = REML.  The best
 |                    default for small-to-moderate N with UN or other
 |                    rich structures -- which is the LaBGAS common case.
 |       DDFM = KR2   Kenward-Roger with the 2009 bias adjustment.
 |       DDFM = SATTERTH   df adjustment only, without the covariance
 |                    correction.  Acceptable; usually close to KR, and
 |                    less prone to convergence trouble.
 |       EMPIRICAL = MBN | MOREL   worth considering INSTEAD when the
 |                    covariance structure itself is doubtful.  KR is
 |                    model-based: it gives a precise df for the model
 |                    you specified, correct or not.
 |       CONTAIN / BETWITHIN / RESIDUAL   fine in large balanced designs,
 |                    where all methods converge.  In an unbalanced
 |                    between-within design, check the df before trusting
 |                    them.
 |
 |     Sanity check on any Type 3 table from a repeated-measures model:
 |     a WITHIN-subject effect should normally carry substantially MORE
 |     denominator df than a BETWEEN-subject effect.  If every effect in
 |     the model shows the same DenDF, something is wrong -- either the
 |     df method is not stratifying, or the df were transcribed.
 |
 |     The macro prints a WARNING when every DenDF is a whole number,
 |     which usually means neither KR nor Satterthwaite was requested.
 |
 |     Kenward MG, Roger JH. Small sample inference for fixed effects
 |     from restricted maximum likelihood. Biometrics 1997;53(3):983-997.
 |     Schaalje GB, McBride JB, Fellingham GW. Adequacy of approximations
 |     to distributions of test statistics in complex mixed linear
 |     models. J Agric Biol Environ Stat 2002;7(4):512-524.
 |     (Verify page numbers before citing.)
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
                        print     = Y,
                        eta2_method = FROM_F,
                        data      = ,
                        class     = ,
                        fixed     = );

    %local _abort _havess _haveprobf _cilevel _mdl _rt _src
           _e2 _o2 _e2lab _o2lab _rtword _srcword _e2m _glmok;
    %let _abort     = 0;
    %let _glmok     = 0;
    %let _e2m       = %upcase(%superq(eta2_method));
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
                    %put WARNING- The two fits do not have the same complete-case set; eta-squared and the F test;
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

                /* How far is DenDF from the observation count?  See
                   ASSUMPTIONS 2: eta-squared from F is inflated by
                   roughly this factor.  Partial eta-squared is immune. */
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
        %if &_havess = 1 %then %do;
            &_e2             = "&_e2lab"
            &_o2             = "&_o2lab"
            mse              = 'Implied MSE (SS_error / DenDF)'
            ss_effect        = 'SS effect (approximate)'
            ss_error         = "SS error (USS of &_rtword residuals)"
            ss_total         = 'SS total (CSS of dependent variable)'
            n_obs            = 'Observations used'
            nobs_dendf       = 'N_obs / DenDF (>1.5 : eta-squared from F unreliable)'
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
                &_o2 = (ss_effect - _df_eff * mse) / (ss_total + mse);
            end;
            else do;
                call missing(ss_effect, &_e2, &_o2);
                eta2_source = 'UNMATCHED';
                put "WARNING: (mixed_effectsize) No Type III SS matched effect " Effect
                    "-- eta-squared suppressed for that row.";
            end;
            label mse       = 'MSE (PROC GLM Type III error)'
                  ss_effect = 'SS effect (PROC GLM Type III)'
                  ss_total  = 'SS total (PROC GLM corrected total)';
            drop _ss_eff _df_eff;
        run;

        %put NOTE: (mixed_effectsize) Eta-squared and omega-squared taken from PROC GLM Type III sums of squares.;
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
    %local _nbad _worst;
    %if &_havess = 1 and &_glmok = 0 %then %do;
        %let _nbad = 0; %let _worst = .;
        proc sql noprint;
            select sum(nobs_dendf > 1.5), max(nobs_dendf)
              into :_nbad trimmed, :_worst trimmed
              from &out where nobs_dendf is not missing;
        quit;
        %if %superq(_nbad) ne and %superq(_nbad) ne . %then %do;
        %if &_nbad > 0 %then %do;
            %put WARNING: (mixed_effectsize) &_nbad effect(s) have N_obs / DenDF > 1.5 (worst: &_worst).;
            %put WARNING- For those rows eta-squared and omega-squared are inflated by roughly that factor,;
            %put WARNING- because MSE is formed as SS_error / DenDF over all N_obs residuals. This is the;
            %put WARNING- expected pattern for a BETWEEN-subject effect in a repeated-measures design.;
            %put WARNING- Report PARTIAL_ETA2 only, or rerun with ETA2_METHOD=DIRECT.;
            %put WARNING- PARTIAL_ETA2 is NOT affected: the MSE cancels out of it.;
        %end;
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
            %put WARNING: (mixed_effectsize) Every DenDF is a whole number.;
            %put WARNING- That usually means the denominator df were left at the default, i.e. neither;
            %put WARNING- DDFM=KR nor DDFM=SATTERTH was requested. The default containment / between-within;
            %put WARNING- methods do not allow for the covariance parameters having been ESTIMATED, and can;
            %put WARNING- be anticonservative in small, unbalanced designs with a rich covariance structure.;
            %put WARNING- Everything in this table depends on DenDF. This is not a claim that KR is always;
            %put WARNING- best -- see ASSUMPTIONS 5 in the macro header for when it is and is not.;
        %end;
    %end;

    /* Same DenDF on every row is a red flag in a between-within design */
    %local _ndistinct _neff;
    proc sql noprint;
        select count(distinct DenDF), count(*) into :_ndistinct trimmed, :_neff trimmed from &out;
    quit;
    %if %superq(_ndistinct) ne and %superq(_neff) ne %then %do;
        %if &_ndistinct = 1 and &_neff > 1 %then %do;
            %put WARNING: (mixed_effectsize) All &_neff effects share the same DenDF.;
            %put WARNING- In a repeated-measures model a WITHIN-subject effect should normally carry;
            %put WARNING- substantially MORE denominator df than a BETWEEN-subject effect. One value for;
            %put WARNING- every effect suggests the df method is not stratifying, or that the df were;
            %put WARNING- transcribed. Check the Type 3 table before reporting anything from here.;
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
                %if &_havess = 1 %then &_e2 &_o2 nobs_dendf;
                eta2_source
                ;
            format FValue 8.2 NumDF 8. DenDF 8.1
                   %if &_haveprobf = 1 %then %str(ProbF pvalue6.4);
                   partial_eta2 partial_eta2_lcl partial_eta2_ucl
                   cohen_f2 8.3
                   %if &_havess = 1 %then %str(&_e2 &_o2 8.4 nobs_dendf 8.2);
                   ;
        run;
        title;
    %end;

    proc datasets library = work nolist nowarn;
        delete _es_ss _es_glm _es_glmoa _es_glmss;
    quit;

    %put NOTE: (mixed_effectsize) Effect sizes written to &out..;
    %put NOTE- Report PARTIAL_ETA2. If you also report &_e2, label it as the output labels do.;

%mend mixed_effectsize;
