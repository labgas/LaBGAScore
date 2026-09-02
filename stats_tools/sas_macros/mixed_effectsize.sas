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
 |  The method previously used in this lab comes from a published SAS
 |  Users Group paper:
 |
 |      Tippey KG, Longnecker MT.  An Ad Hoc Method for Computing
 |      Pseudo-Effect Size for Mixed Models.  Texas A&M University.
 |
 |  Its macro computes three statistics -- eta_2, omega_2 and
 |  partial_eta_2 -- into one dataset, so it is easy to report one of them
 |  under the name of another.  In the ME/CFS analysis the reported
 |  quantity turned out to be partial_eta_2 itself: the NAME was right and
 |  the arithmetic was wrong, which is the harder failure to spot.  It
 |  forms the error sum of squares as
 |
 |      ss_error = mse * (_freq_ - numdf);
 |
 |  which uses the NUMBER OF OBSERVATIONS minus the effect's NUMERATOR
 |  degrees of freedom, where the model's DENOMINATOR degrees of freedom
 |  belong.  That line is verbatim in the published macro's appendix, so
 |  this is a defect in the method as published, not a local slip -- and
 |  any analysis built on that paper inherits it.
 |
 |  WHEN IT BITES, AND WHEN IT DOES NOT.  The mse cancels, leaving
 |
 |      partial_eta_2 = qF / (qF + _freq_ - q)
 |
 |  against the correct qF / (qF + DenDF).  The two AGREE if and only if
 |  _freq_ - q = DenDF, and that needs TWO independent things to hold:
 |  the residuals dataset must contain only the rows the model fitted,
 |  and DenDF must be the observation-scale residual df.
 |
 |  The second is decided by the DF METHOD in force -- NOT by whether the
 |  effect crosses subjects:
 |
 |    CONTAINMENT, falling through -- a RANDOM statement is present and
 |      no random effect contains the fixed effect
 |        -> N - rank(XZ) for EVERY effect.  AGREES.
 |    BETWEEN-WITHIN with a residual variance term (CS, VC, AR(1))
 |        -> stratifies.  Agrees for within-subject effects, FAILS for
 |           between-subject ones.
 |    BETWEEN-WITHIN under type = UN
 |        -> participant scale for EVERY effect, because a saturated
 |           R-side leaves no within-subject stratum.  FAILS for all,
 |           within-subject effects included.
 |    KR / KR2 / SATTERTH
 |        -> fractional, in general neither scale.  FAILS.
 |
 |  The published paper's worked example is the SECOND case, which is why
 |  the between/within distinction looks like the governing one.  Its
 |  Case 1 (fev1: 72 patients x 8 hours = 576 observations) reports
 |  partial eta-squared .029 / .082 / .057 for Drug / Hour / Drug*Hour;
 |  recomputed against the design's real denominator df, the two
 |  WITHIN-subject effects move by 1.2x and the BETWEEN-subject Drug
 |  effect by 6.9x (.029 -> .199).  Do not generalise from it: a project
 |  whose models all carry a RANDOM statement lands in the FIRST case and
 |  agrees to 1-2% even for between-subject factors.
 |
 |  So when auditing an old table the question is not "was this script
 |  used", nor "does the effect cross subjects", but "which DF METHOD was
 |  in force, and did the residuals dataset hold only fitted rows".
 |
 |  ---------------------------------------------------------------------
 |  OBSERVATIONS READ versus OBSERVATIONS USED
 |  ---------------------------------------------------------------------
 |  An OUTPM= dataset has one row per observation the procedure READ, not
 |  per observation it USED.  Rows dropped for a missing response or
 |  covariate are still there, carrying a missing residual.  In a long
 |  file laid out on the longest outcome's timepoint grid the gap can be
 |  enormous: in the analysis this macro replaces, PROC MIXED read 1384
 |  rows for every model while using between 483 and 1312 of them, and
 |  the old code's _FREQ_ took the value 1384 in all of them.
 |
 |  This macro is immune to that by construction:
 |    * PARTIAL_ETA2, PARTIAL_OMEGA2, PARTIAL_EPSILON2 and the confidence
 |      interval read ONLY the Type 3 table.  DenDF already reflects the
 |      observations the model actually fitted.  Nothing to get wrong.
 |    * ETA2 and OMEGA2 filter the residuals dataset to rows with a
 |      non-missing residual before taking any sum of squares, so
 |      SS_total, SS_error and N_OBS all describe the same set of rows --
 |      the fitted ones.  N_OBS is therefore observations USED.
 |  The macro prints a NOTE giving both counts whenever they differ, so
 |  the size of the gap is visible rather than silently absorbed.
 |
 |  ---------------------------------------------------------------------
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
 |  --- degrees-of-freedom diagnostics (all optional) -------------------
 |  modelinfo  Dataset from ODS OUTPUT MODELINFO=.  Supplies the
 |             covariance structure, the df method actually used, the
 |             residual variance method and the subject effect, so the
 |             macro can tailor its df checks instead of guessing.
 |  dimensions Dataset from ODS OUTPUT DIMENSIONS=.  Supplies Columns in
 |             Z (whether a RANDOM statement is active) and the subject
 |             count -- but see NSUBJECTS=.
 |  nsubjects  Number of subjects.  Give this explicitly whenever you
 |             want the between-subject df check.  It is needed because
 |             the Dimensions subject count is NOT always usable: a
 |             RANDOM effect with no SUBJECT= option spans subjects, so
 |             SAS reports Subjects=1 and Max Obs per Subject = N.  The
 |             macro detects that case and asks for this parameter.
 |  between    Space-separated list of effect names that are
 |             BETWEEN-subject, exactly as they appear in TESTS3.  With
 |             NSUBJECTS= this turns the df check from a hint into a
 |             statement: a between-subject effect cannot have more
 |             denominator df than there are subjects.
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
 |  2. *** ETA-SQUARED RECONSTRUCTED FROM F IS NOT TRUSTWORTHY FOR A
 |     MIXED OR MARGINAL MODEL.  USE ETA2_METHOD = DIRECT. ***
 |
 |     SS_effect = NumDF x FValue x MSE is exact in a fixed-effects
 |     ANOVA, where F really is MS_effect / MS_error and the error df
 |     really is N - rank(X).  Neither holds here.  Two separate things
 |     go wrong:
 |
 |       (a) MSE is formed as SS_error / DenDF, but SS_error is the USS
 |           of the residuals over ALL N_obs observations.  When DenDF
 |           is participant-level -- a BETWEEN-subject effect in a
 |           repeated design -- that MSE is inflated by about
 |           df_resid / DenDF.
 |       (b) F from PROC MIXED is a Wald statistic on GLS estimates,
 |           not a ratio of orthogonal sums of squares.  It can be much
 |           larger or much smaller than the corresponding OLS F.
 |
 |     The net error is the PRODUCT of the two, and they can compound or
 |     partly cancel.  Verified against R (effectsize, lme4/lmerTest) and
 |     Python (pingouin, statsmodels) on a simulated 120-subject x
 |     4-timepoint design:
 |
 |       effect        N_obs/DenDF   MSE inflation   actual eta2 error
 |       group             4.10          4.03             1.29x
 |       time              1.36          1.34             4.47x
 |       group x time      1.36          1.34             4.47x
 |
 |     So the error is NOT predictable from the degrees of freedom, and
 |     N_OBS/DenDF is a diagnostic of (a) only, not a correction factor.
 |     An earlier version of this file claimed the inflation was roughly
 |     N_obs/DenDF.  That was wrong; the table above is the evidence.
 |
 |     NOBS_DENDF is still reported because it makes (a) visible, but the
 |     macro now warns on EVERY row whenever ETA2_METHOD = FROM_F, since
 |     the problem is not confined to rows with a small DenDF.
 |
 |     PARTIAL_ETA2, PARTIAL_OMEGA2 and PARTIAL_EPSILON2 ARE NOT
 |     AFFECTED -- they use only NumDF, DenDF and FValue.  This caveat is
 |     about eta-squared and (classical) omega-squared alone.
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
 |  5. Choose the denominator-df method deliberately.  Everything here
 |     depends on DenDF, so this is the single most consequential choice
 |     upstream of the macro.
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
 |     Why a fixed standard: otherwise SAS picks one of two defaults FOR
 |     you, decided by whether a RANDOM statement happens to be present,
 |     and the two are on different SCALES --
 |
 |       RANDOM present  -> CONTAIN    falls through to N - rank(XZ),
 |                                     i.e. OBSERVATION-scale df
 |       REPEATED only   -> BETWITHIN  SUBJECT-scale df for a
 |                                     between-subject predictor
 |
 |     so the same data and the same F can give a partial eta-squared
 |     several times smaller in one model than in the other.  Naming the
 |     method removes that accident and buys the small-sample covariance
 |     correction as well.  Full reasoning, the SAS quotes and a worked
 |     contrast: README.md, "Choose a denominator-df method".
 |
 |     THE OPTIONS, AND HOW TO WRITE THEM
 |     ................................................................
 |     ddfm = kr2   Kenward-Roger with the Kenward & Roger (2009)
 |                  precision estimator.  THE STANDARD.  Does two
 |                  things: inflates the covariance matrix of the fixed
 |                  effects to allow for the covariance parameters
 |                  having been ESTIMATED, and derives a Satterthwaite-
 |                  type df from it.  So it changes both F and DenDF.
 |     ddfm = kr    Kenward & Roger (1997).  IDENTICAL to kr2 for a
 |                  covariance structure that is linear in its
 |                  parameters -- UN, CS, VC, TOEP, the LaBGAS common
 |                  case.  The two differ only for NONLINEAR structures
 |                  (AR(1), ARH(1), CSH, TOEPH, SP()), where the KR
 |                  second-order term can SHRINK standard errors and is
 |                  not invariant to reparameterisation.  Use kr only
 |                  where kr2 is unavailable.
 |     ddfm = kr(firstorder)   drops the second-derivative term.
 |                  Subsumed by kr2; no reason to choose it.
 |     ddfm = satterth   df adjustment only, no covariance correction.
 |                  Acceptable fallback: usually close to KR, and less
 |                  prone to convergence trouble.
 |     ddfm = contain | betwithin | residual   the defaults.  Fine in
 |                  large balanced designs, and fine under type=UN (see
 |                  below).  Not wrong -- but they do not allow for the
 |                  covariance parameters having been estimated, and
 |                  they leave the choice to SAS.
 |     empirical = mbn | morel   goes on the PROC MIXED statement, not
 |                  the MODEL statement.  Consider it INSTEAD of KR when
 |                  the covariance STRUCTURE itself is doubtful: KR is
 |                  model-based, and gives a precise df for the model
 |                  you specified, right or wrong.
 |
 |     Three conditions on KR / KR2:
 |       * they are defined for METHOD = REML (the PROC MIXED default).
 |         A likelihood-ratio test on FIXED effects needs ML, where KR
 |         does not apply -- do not mix the two in one table.
 |       * if a variance component is estimated at the boundary (0 in
 |         Covariance Parameter Estimates), the adjustment is built on a
 |         Hessian at that boundary.  Check before trusting the df.
 |       * a RANDOM effect with no SUBJECT= spans subjects, so SAS
 |         cannot block the problem and Dimensions reports Subjects = 1.
 |         Fix the blocking before reaching for KR -- and note the macro
 |         cannot use that subject count, which is why NSUBJECTS=
 |         exists.
 |
 |     UNDER type = UN, DO NOT "FIX" THE DEFAULT.  BETWITHIN prints the
 |     same DenDF for every effect there, and that is very nearly right:
 |     the exact multivariate df are n - rank for between-subject
 |     effects and n - rank - q + 1 for within-subject ones, so SAS is
 |     out by a couple of df (161 against 159 in a verified 164 x 4
 |     design; at most 0.006 on partial eta-squared).  KR2 still earns
 |     its place there -- for the covariance correction, not the df --
 |     and will move F, and so the effect sizes, slightly DOWNWARD.
 |     Do NOT expect a within-subject df of n*t - n - rank; that is the
 |     compound-symmetry answer and does not apply under UN.
 |
 |     THE MACRO RUNS FINE ON ANY OF THESE, defaults included.  Nothing
 |     about the df method can stop it -- the df diagnostics are NOTE
 |     and WARNING text only.  What changes is the INTERPRETATION:
 |     partial eta-squared is conditioned on DenDF and has no df-free
 |     value, so report the df method beside the effect size.
 |
 |     Kenward MG, Roger JH. Small sample inference for fixed effects
 |     from restricted maximum likelihood. Biometrics 1997;53:983-997.
 |     Kenward MG, Roger JH. An improved approximation to the precision
 |     of fixed effects from restricted maximum likelihood.
 |     Comput Stat Data Anal 2009;53:2583-2595.
 |     Schaalje GB, McBride JB, Fellingham GW. Adequacy of
 |     approximations to distributions of test statistics in complex
 |     mixed linear models. J Agric Biol Environ Stat 2002;7:512-524.
 |     (Verify page numbers before citing.)
 |
 |  ---------------------------------------------------------------------
 |  CROSS-VALIDATED AGAINST R AND PYTHON  (2026-09-01)
 |  ---------------------------------------------------------------------
 |  Simulated repeated-measures design, 120 subjects x 4 timepoints,
 |  between-subject GROUP, within-subject TIME, subject-level covariate,
 |  random intercept.  Fitted with lme4/lmerTest (Kenward-Roger) and with
 |  pingouin's mixed_anova.
 |
 |  PARTIAL ETA-SQUARED -- exact agreement, three implementations:
 |      this macro          (NumDF x F) / (NumDF x F + DenDF)
 |      R  effectsize::F_to_eta2()          identical to 1e-16
 |      Py pingouin.mixed_anova(effsize='np2')  identical to 1e-16
 |  pingouin computes it from a REAL repeated-measures sum-of-squares
 |  decomposition, not from F, and still lands on the same number.  The
 |  formula is an algebraic identity, not an approximation.
 |
 |  CONFIDENCE INTERVAL -- this macro is the more accurate of the two.
 |  Both invert the noncentral F.  Checking where the returned bounds
 |  actually put the tail probability (target 0.975 / 0.025):
 |      this macro (FNONCT, exact inversion)  0.975000 / 0.025000
 |      R effectsize (optim-based NCP search) 0.972129 / 0.020694
 |  Bounds differ by about 0.003-0.004 in eta-squared units.  DO NOT
 |  "correct" this macro to reproduce R's numbers.
 |
 |  ETA-SQUARED -- see ASSUMPTIONS 2.  Neither FROM_F variant reproduces
 |  a real Type III decomposition; ETA2_METHOD = DIRECT does, by
 |  construction, and was checked against car::Anova(type=3) in R and
 |  statsmodels.anova_lm(typ=3) in Python (both with sum-to-zero
 |  contrasts, which Type III requires).
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
                        fixed     = ,
                        modelinfo = ,
                        dimensions = ,
                        nsubjects = ,
                        between   = );

    %local _abort _havess _haveprobf _cilevel _mdl _rt _src
           _e2 _o2 _e2lab _o2lab _rtword _srcword _e2m _glmok
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
        %let _o2lab   = %str(Omega-squared, CLASSICAL/non-partial, proportion of total variance (approximate) -- NOT the same as PARTIAL_OMEGA2);
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
            %put NOTE: (mixed_effectsize) &resids holds &_n_rows rows; &_n_obs have a residual.;
            %put NOTE- Sums of squares and N_OBS are taken over the &_n_obs rows the model actually;
            %put NOTE- used. The other rows were read but not fitted (missing response or covariate,;
            %put NOTE- or timepoint slots belonging to a different outcome in a long file).;
        %end;

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

            /* Partial omega-squared and partial epsilon-squared. Like
               partial eta-squared these need ONLY NumDF, DenDF and
               FValue, so a reader can recompute them from the same
               table -- and both are less positively biased than partial
               eta-squared. Definitions match R's effectsize package
               (F_to_omega2, F_to_epsilon2), verified numerically.
               They can go negative when F < 1; that is expected, and
               such a value is conventionally read as 0.            */
            partial_omega2   = (NumDF * (FValue - 1)) / (NumDF * FValue + DenDF + 1);
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
            partial_omega2   = 'Partial omega-squared (less biased; from F and df only)'
            partial_epsilon2 = 'Partial epsilon-squared (less biased; from F and df only)'
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
        %put WARNING- is NOT predictable from the df. Report PARTIAL_ETA2 / PARTIAL_OMEGA2 /;
        %put WARNING- PARTIAL_EPSILON2, which are unaffected, or rerun with ETA2_METHOD=DIRECT.;
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
        title3 "&_e2 computed from &_rtword residuals -- an approximation for a mixed model; see the macro header";
        %end;

        proc print data = &out noobs label;
            var Effect NumDF DenDF FValue
                %if &_haveprobf = 1 %then ProbF;
                partial_eta2 partial_eta2_lcl partial_eta2_ucl
                partial_omega2 partial_epsilon2 cohen_f2
                %if &_havess = 1 %then &_e2 &_o2 nobs_dendf;
                %if &_havesub = 1 %then dendf_ratio;
                %if %superq(between) ne %then is_between;
                eta2_source
                ;
            format FValue 8.2 NumDF 8. DenDF 8.1
                   %if &_haveprobf = 1 %then %str(ProbF pvalue6.4);
                   partial_eta2 partial_eta2_lcl partial_eta2_ucl
                   partial_omega2 partial_epsilon2 cohen_f2 8.3
                   %if &_havess = 1 %then %str(&_e2 &_o2 8.4 nobs_dendf 8.2);
                   %if &_havesub = 1 %then %str(dendf_ratio 8.2);
                   ;
        run;
        title;
    %end;

    proc datasets library = work nolist nowarn;
        delete _es_ss _es_glm _es_glmoa _es_glmss;
    quit;

    %put NOTE: (mixed_effectsize) Effect sizes written to &out..;
    %put NOTE- Report PARTIAL_ETA2 (optionally with PARTIAL_OMEGA2 / PARTIAL_EPSILON2, which are less biased);
    %put NOTE- and are recomputable by a reader from the F and df in the same table.;
    %put NOTE- If you also report &_e2, label it as the output labels do; it is the CLASSICAL, non-partial;
    %put NOTE- omega/eta family and is not comparable with PARTIAL_OMEGA2.;

%mend mixed_effectsize;
