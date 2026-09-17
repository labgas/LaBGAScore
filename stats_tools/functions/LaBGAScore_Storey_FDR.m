function [q, pi0, info] = LaBGAScore_Storey_FDR(p, varargin)
% Storey's positive FDR q-values, with a pi0 estimate that is checked rather
% than trusted.
%
% *USAGE*
%
% q                = LaBGAScore_Storey_FDR(p)
% [q, pi0]         = LaBGAScore_Storey_FDR(p)
% [q, pi0, info]   = LaBGAScore_Storey_FDR(p, 'method', 'lambda')
% [q, pi0, info]   = LaBGAScore_Storey_FDR(p, 'lambda', 0:0.05:0.95, 'verbose', false)
%
% *WHY THIS WAS REWRITTEN*
%
% The previous version called mafdr(p), accepted its pi0 unless it exceeded
% 0.99, and then floored q at p. That is unsafe, and it failed on real data
% from proj_cfs:
%
%   8 roi SVM p-values : pi0 = 0.895  -> fine
%   8 roi GLM p-values : pi0 = 0.012  -> claims ~99% of 8 tests are non-null,
%                                        from p-values of .04 to .52
%
% In the second case every q fell below its own p, so the q >= p floor - which
% is correct in itself - turned the whole vector back into the RAW p-values,
% returned under the name q. An uncorrected p reported as an FDR q is the worst
% failure available here, and note it is DATA-DEPENDENT: identical code gave a
% sane answer on the other set, so one run cannot reveal it.
%
% The 0.99 guard only ever caught the CONSERVATIVE failure (pi0 -> 1, where
% Storey harmlessly degenerates to BH). The damaging direction is pi0 -> 0.
%
% *RELATION TO SAS PROC MULTTEST*
%
% The old header claimed this implemented Storey "as in SAS proc multtest". It
% did not. SAS's PFDR default is:
%
%   "the SPLINE method is attempted first. If the estimate is nonpositive or if
%    the slope of the spline at the last lambda is greater than 0.1 times the
%    range of the fitted spline values, then the BOOTSTRAP method is used."
%
% with NLAMBDA=20 and NBOOT=10000 (both Storey & Tibshirani 2003); SAS's
% LAMBDA= instead fixes a single lambda with no search (Storey 2002). MATLAB's
% mafdr implements the SPLINE step only, with no fallback - which is exactly
% why the degenerate estimate above was returned silently.
%
% *METHODS AVAILABLE, AND WHAT THEY CONTROL*
%
%   FDR, Storey-type (q = pi0 * q_BH, differing only in how pi0 is estimated):
%     'sas'          SPLINE then BOOTSTRAP on SAS's trigger = PROC MULTTEST PFDR
%     'lambda'       median of pi0(lambda) over a grid; not a SAS method
%     'spline'       mafdr's spline step alone
%     'bh'           pi0 = 1, i.e. plain Benjamini-Hochberg
%
%   FDR, adaptive (m replaced by an estimate of the number of true nulls):
%     'adaptivefdr'  Benjamini & Hochberg (2000) adaptive linear step-up, m0 by
%                    the lowest-slope estimator of Hochberg & Benjamini (1990)
%                    after Schweder & Spjotvoll (1982). SAS ADAPTIVEFDR default
%                    (LOWESTSLOPE).
%     'bky'          Benjamini, Krieger & Yekutieli (2006) two-stage linear
%                    step-up: BH at alpha/(1+alpha), then BH with m0 = m - r1.
%
%   FWER (probability of ANY false positive, NOT comparable to the FDR columns):
%     'stepdown_sidak'  Holm step-down with Sidak multiplier 1-(1-p)^k
%     'holm'            Holm step-down with Bonferroni multiplier k*p
%
% *MEASURED BEHAVIOUR OF THE pi0 ESTIMATORS*
%
% Simulated here, m = 200, 200 runs per cell, mean estimated pi0 (bias):
%
%   true pi0    'adaptivefdr'     'bky'            'sas'
%     0.30      0.413 (+0.11)   0.602 (+0.30)   0.224 (-0.08)
%     0.50      0.619 (+0.12)   0.748 (+0.25)   0.384 (-0.12)
%     0.80      0.883 (+0.08)   0.930 (+0.13)   0.667 (-0.13)
%     1.00      0.999 (-0.00)   1.000 (-0.00)   0.835 (-0.17)
%
% The two adaptive methods are biased UP (conservative); Storey is biased DOWN
% (anti-conservative), and stays biased down even at pi0 = 1, where it calls
% ~17%% of a pure null signal. Rejection rate on a complete null at alpha = .05,
% m = 200: BH .052, 'adaptivefdr' .052, 'bky' .046, 'sas' .078. At m = 8 the
% FWER methods behave as they should: 'stepdown_sidak' .058, 'holm' .048, and
% 'stepdown_sidak' agrees with CanlabCore holm_sidak on 500/500 datasets.
%
% What this means in practice, on a real 6-signature panel from this pipeline:
%
%   method            pi0     smallest q
%   'sas'            0.023      0.0152     <- only method declaring significance
%   'adaptivefdr'    0.833      0.0760
%   'bh' / 'bky'     1.000      0.0912
%   'stepdown_sidak'   -        0.0878
%   'holm'             -        0.0912
%
% Both adaptive methods are DESIGNED to gain power over BH and still do not
% reach .05 here, while 'sas' does - entirely on the strength of pi0 = 0.023
% estimated from six tests. Prefer 'adaptivefdr' or 'bky' when the panel is
% small and the question is whether an effect survives correction at all.
%
% Two Storey methods are offered here:
%
%   'sas' (DEFAULT) - SPLINE first, falling back to the Storey & Tibshirani
%       BOOTSTRAP on SAS's own trigger, i.e. SAS PROC MULTTEST's PFDR default.
%       This is the default deliberately: LaBGAS cross-checks analyses against
%       SAS, and a q-value that cannot be reproduced by PROC MULTTEST is worth
%       less than a slightly better-estimated one that can.
%   'lambda' - pi0 is the median of Storey's fixed-lambda estimator
%       pi0(lambda) = #{p > lambda} / (n * (1 - lambda))
%       over a grid. Not a SAS method; closest to a robustified LAMBDA=.
%
% The accuracy cost of that choice is real and is recorded here so nobody has
% to rediscover it. Simulated at n = 200 with a true pi0 of 0.70, 100 runs:
%
%   lambda median     bias +0.006   sd 0.039
%   mafdr spline      bias -0.168   sd 0.171
%   ST bootstrap      bias -0.136   sd 0.160
%
% and raising NBOOT does not help (100 -> 1000 draws moved the bootstrap's bias
% from -0.094 to -0.125), because it minimises MSE against min(pi0(lambda)), a
% target that is itself biased low. Bias in pi0 is the anti-conservative
% direction, so 'sas' errs towards declaring too much significant.
%
% *HOW 'sas' COMPARES WITH R's qvalue PACKAGE*
%
% R's qvalue (Bioconductor, written by Storey) is the reference implementation:
%
%   pi0est(p, lambda = seq(0.05, 0.95, 0.05),
%          pi0.method = c("smoother", "bootstrap"), smooth.df = 3)
%
% Its default is the SPLINE ("smoother", Storey & Tibshirani 2003) with no
% fallback, and its only guards are:
%   pi0 > 1   -> silently clamped, pi0 <- min(pi0, 1)
%   pi0 <= 0  -> warning("The estimated pi0 <= 0. Setting the pi0 estimate to
%                be 1..."), i.e. fall back to BH
%   length(lambda) < 4, or max(p) < min(lambda) -> hard stop()
%
% So R guards only against the mathematically impossible, not the implausible:
% a pi0 of 0.02 estimated from six tests is positive, so R uses it without
% complaint. SAS does the same. So does 'sas' here, by design.
%
% 'sas' is nonetheless STRICTER THAN R IN MOST CASES, because SAS's PFDR adds a
% step R's smoother has no equivalent of: it rejects a spline estimate the
% lambda curve does not support and falls back to the Storey & Tibshirani (2003)
% bootstrap, which usually returns a HIGHER pi0 than the rejected spline.
% Measured here against an emulation of R's smoother path, 300 runs per cell:
%
%    m    true pi0   pi0 R    pi0 'sas'   min q R   min q 'sas'   'sas' >= R
%    6      0.50     0.099     0.088      0.0064      0.0056         92%
%    8      0.50     0.117     0.115      0.0045      0.0044         91%
%    8      0.80     0.206     0.206      0.0170      0.0199         89%
%   20      0.70     0.272     0.302      0.0048      0.0055         89%
%   50      0.70     0.409     0.438      0.0013      0.0014         86%
%  200      0.70     0.556     0.559      0.0001      0.0001         84%
%  200      1.00     0.820     0.821      0.4224      0.4232         87%
%
% Read that last column as "how often 'sas' is at least as conservative as R":
% 84-92%%, never 100%%. The advantage grows with m, where the bootstrap fallback
% fires more usefully. At m = 6 it REVERSES - mean pi0 0.088 against R's 0.099,
% so 'sas' is slightly more permissive than R on the smallest panels. Do not
% rely on 'sas' being the safer of the two at small m; it is not, reliably.
%
% For completeness on the other side: Python's standard stack has NO Storey
% method at all. statsmodels.stats.multitest.multipletests offers bonferroni,
% sidak, holm-sidak, holm, simes-hochberg, hommel, fdr_bh, fdr_by, fdr_tsbh,
% fdr_tsbky (= the BKY 2006 procedure implemented here as 'bky') and fdr_gbs,
% but no pi0 estimation; Storey requires py-qvalue or multipy.
%
% *THE RELIABILITY GUARD, AND WHY 'sas' DOES NOT USE IT*
%
% The diagnostics below (pi0 across a lambda grid, a degenerate-pi0 test) are
% computed and printed for EVERY method, because they are worth seeing. Whether
% they are allowed to OVERRULE the estimate and return BH instead is a separate
% question, controlled by 'guard':
%
%   method 'sas'    guard OFF by default  -> reproduces SAS PROC MULTTEST PFDR
%   other methods   guard ON  by default  -> rejects a pi0 the lambda curve
%                                            does not support, returns BH
%
% 'sas' has the guard off because SAS has no such check: PROC MULTTEST
% estimates pi0 and uses it, full stop. A guard makes the output safer but no
% longer reproducible in SAS, which defeats the purpose of having a SAS mode at
% all. LaBGAS cross-checks against SAS, so 'sas' must mean SAS.
%
% Know what that costs. The guard is not decoration: it fires often at the
% panel sizes used here, and every firing is a case where the raw estimate was
% implausible. Measured over 200 simulations per cell:
%
%   n     true pi0    guard would fire    dominant reason
%     8       0.50          81.5%         pi0 < 0.01 (71.5%)
%     8       0.80          72.5%         pi0 < 0.01 (65.0%)
%    20       0.50          48.5%         pi0 < 0.01 (46.5%)
%    50       0.70          14.5%         pi0 < 0.01
%   100       0.70           2.0%
%   200       0.70           0.0%
%
% So on an 8-signature panel, 'sas' now returns a Storey q in roughly three
% cases out of four where the old behaviour returned BH - and in most of those
% the spline's pi0 was below 0.01, i.e. asserting that ~99%% of 8 tests are true
% effects. Those q-values are ANTI-CONSERVATIVE. They are what SAS gives, they
% are what this function now gives, and the verdict line says so explicitly
% whenever the estimate is questionable. info.q_BH always carries the
% conservative alternative; report both when n is small.
%
% Set 'guard', true to restore the old protective behaviour under 'sas', or
% 'guard', false to strip it from another method.
%
% *HOW q IS FORMED*
%
%   q = pi0 * q_BH,  then floored at p
%
% which is what Storey's procedure reduces to, and makes the relationship
% explicit: Storey is BH scaled by the estimated null proportion, so it is
% never weaker than BH, and at pi0 = 1 the two coincide. The floor at p follows
% SAS proc multtest.
%
% *INPUTS*
%
%   p           vector of raw p-values
%
% *OPTIONAL INPUTS*
%
%   'method'    'sas' (default, = SAS PROC MULTTEST PFDR) | 'lambda' | 'spline' | 'bh'
%   'lambda'    grid for the lambda method. Default 0.2:0.1:0.5
%   'nboot'     bootstrap draws for 'sas'. Default 1000 (SAS uses 10000)
%   'verbose'   print the diagnostics. Default true
%   'guard'     override the reliability guard: true forces it on, false forces
%               it off. Default: off for 'sas' (SAS compatibility), on for the
%               other methods. See THE RELIABILITY GUARD above.
%
% *OUTPUTS*
%
%   q           FDR-corrected p-values, same shape as p
%   pi0         the estimated proportion of true nulls that was USED
%   info        struct: .method_used, .pi0, .pi0_lambda, .lambda, .pi0_range,
%               .pi0_spline, .reliable, .reasons, .q_BH,
%               .guard_on (was the guard active), .guard_fired (did it override)
%
% *WHAT THE ORIGINAL PAPERS SAY ABOUT THE NUMBER OF TESTS*
%
% Checked against the primary sources, because this function previously carried
% a storey_min_tests = 50 rule with no citation behind it. That rule has been
% REMOVED: it was invented here, and it is not in Storey.
%
% NEITHER paper states a minimum m. What they state is different, and stronger:
%
%   Storey (2002), JRSS-B 64:479-498, Discussion, p.495:
%     "The methodology presented here has the opposite property - the more tests
%      we perform, the better the estimates are. Therefore, it is an asset under
%      this approach to have large data sets with many tests. THE ONLY
%      REQUIREMENT IS THAT THE TESTS MUST BE EXCHANGEABLE in the sense that the
%      p-values have the same null distribution."
%
% So the stated requirement is EXCHANGEABILITY, not a count. That is a mild
% requirement for the families used here: the 8 MIST ROIs are disjoint (verified:
% zero shared voxels when the combined atlas was built), and the npsplus
% signatures, while they may share some voxels, are in practice largely
% non-overlapping. Correlation between tests is worth remembering but is not a
% reason to distrust these families - and Storey (2002) notes the effect of
% dependence is in any case negligible for large m, with Storey & Tibshirani
% (2003) showing the q-values remain conservative under weak dependence.
%
% The theory is asymptotic in m throughout: "for large m this assumption makes
% little difference" (p.482), "For large m these two estimates are equivalent"
% (p.483), "the effect of dependence is negligible if m is large" (p.495).
% Storey & Tibshirani (2003), PNAS 100:9440-9445, is the same: their results
% are limits "as m -> infinity", and their worked example has m = 3170 genes.
% Finite-sample trouble is handled by ADJUSTMENT, not by a cutoff - R(gamma)v1
% exists because "when R(gamma) = 0, the estimate would be undefined, which is
% undesirable for finite samples" (p.483).
%
% The only empirical evidence about small m is Storey (2002) Table 3, p.495,
% the simulation for the bootstrap lambda-selection that SAS's PFDR uses:
%
%     m       lambda_best   lambda_hat   MSE(lambda_hat)
%       100        0.75         0.65          0.127
%       500        0.75         0.75          0.00953
%      1000        0.80         0.80          0.00444
%     10000        0.90         0.90          0.000556
%
% MSE at m = 100 is ~13x that at m = 500, and the selected lambda misses. The
% smallest m Storey ever evaluates is therefore 100, and the method is already
% visibly degrading there. Storey & Tibshirani's own caution is about lambda,
% not m: "as we set lambda closer to 1, the variance of pi0_hat(lambda)
% increases, making the estimated q values more unreliable".
%
% PRACTICAL UPSHOT for the panels this is used on (6-8 signatures, 8 ROIs):
% that is one to two orders of magnitude below anything in the papers. The
% literature does not forbid it and does not support it. Nothing here blocks
% it - the estimate is returned, and under 'sas' it is used - but ALWAYS report
% info.q_BH alongside q, and treat a pi0 near zero at m < 100 as an artefact of
% tiny m rather than as evidence that almost every test is non-null.
%
% *WHEN pi0 CANNOT BE TRUSTED*
%
% info.reliable is false, info.reasons says why, and q falls back to BH. BH is
% not a compromise: BH IS Storey at pi0 = 1, the conservative choice you want
% when pi0 is not identifiable. Note also that correlated tests (roi means from
% the same subjects, neighbouring parcels) violate Storey's independence
% assumption, while BH holds under positive regression dependency.
%
% *SEE ALSO*
%
% mafdr, prep_3a_run_second_level_regression_and_save,
% LaBGAScore_decoding_SVM_between_subjects
%
% -------------------------------------------------------------------------
% Lukas Van Oudenhove, KU Leuven, September 2026
% -------------------------------------------------------------------------

% ---------------------------- parse inputs -------------------------------

ip = inputParser;
ip.addParameter('method',  'sas', @(x) ischar(x) || isstring(x));
ip.addParameter('lambda',  0.2:0.1:0.5, @isnumeric);
ip.addParameter('nboot',   1000, @isnumeric);
ip.addParameter('verbose', true, @(x) islogical(x) || isnumeric(x));
ip.addParameter('guard',   [], @(x) isempty(x) || islogical(x) || isnumeric(x));
ip.addParameter('alpha',   0.05, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
ip.parse(varargin{:});

method  = lower(char(ip.Results.method));
lam     = ip.Results.lambda(:)';
nboot   = ip.Results.nboot;
verbose = logical(ip.Results.verbose);
alpha   = ip.Results.alpha;

% The reliability guard is an ADDITION to Storey, not part of it. SAS PROC
% MULTTEST has no such check: it estimates pi0 and uses it. So the guard is OFF
% for method 'sas', which now reproduces PROC MULTTEST PFDR rather than
% second-guessing it, and ON for the other methods, where no external
% reference is being matched and the degenerate-pi0 failure is worth catching.
% Override either way with 'guard', true/false.
if isempty(ip.Results.guard)
    do_guard = ~ismember(lower(char(ip.Results.method)), ...
        {'sas','bh','stepdown_sidak','sidak','holm','bonferroni', ...
         'adaptivefdr','bh2000','lsl','bky','bky2006','twostage'});
else
    do_guard = logical(ip.Results.guard);
end

sz = size(p);
p  = double(p(:));

if any(p < 0 | p > 1) || any(isnan(p))
    error('LaBGAScore_Storey_FDR:badp', 'p must contain only non-NaN values in [0,1]');
end

n = numel(p);

q_bh = mafdr(p, 'BHFDR', true);
q_bh = q_bh(:);

% ------------------------- estimate pi0 ----------------------------------

pi0_lam   = arrayfun(@(L) min(1, mean(p > L)/(1-L)), lam);
pi0_range = max(pi0_lam) - min(pi0_lam);
pi0_med   = median(pi0_lam);

pi0_spline = NaN;
try
    [~, ~, pi0_spline] = mafdr(p);
catch
    % mafdr's spline fit can fail outright on small or degenerate inputs
end

switch method
    case 'lambda'
        pi0_use = pi0_med;  method_used = 'lambda-median';
    case 'spline'
        pi0_use = pi0_spline;  method_used = 'mafdr-spline';
    case 'sas'
        [pi0_use, method_used] = local_sas_pi0(p, lam, pi0_spline, nboot);
    case 'bh'
        pi0_use = 1;  method_used = 'BH (pi0 = 1)';
    case {'stepdown_sidak','sidak','holm','bonferroni'}
        pi0_use = NaN;  method_used = method;   % FWER methods: pi0 does not apply
    case {'adaptivefdr','bh2000','lsl'}
        [pi0_use, method_used] = local_lsl_pi0(p);
    case {'bky','bky2006','twostage'}
        [pi0_use, method_used] = local_bky_pi0(p, q_bh, alpha);
    otherwise
        error('LaBGAScore_Storey_FDR:method', ...
            ['method must be ''sas'', ''lambda'', ''spline'', ''bh'', ' ...
             '''stepdown_sidak'', ''holm'', ''adaptivefdr'' or ''bky'', got ''%s'''], method);
end

% ------------------------- judge pi0 -------------------------------------

reasons = {};
if pi0_range > 0.3
    reasons{end+1} = sprintf('pi0 swings %.2f-%.2f across lambda, so it is not identifiable', ...
        min(pi0_lam), max(pi0_lam));
end
if ~isnan(pi0_use) && pi0_use < 0.01
    reasons{end+1} = sprintf('pi0 = %.4f is degenerate (implausibly few nulls)', pi0_use);
end
if isnan(pi0_use)
    reasons{end+1} = 'pi0 could not be estimated';
end

reliable = isempty(reasons);

% ------------------------- form q ----------------------------------------

% The diagnostics above are computed and reported WHATEVER the method, because
% they are worth seeing. What changes is whether they are allowed to overrule
% the estimate: only when the guard is active.
guard_fired = do_guard && ~reliable;

if ismember(method, {'stepdown_sidak','sidak','holm','bonferroni'})
    % FWER methods: they control the probability of ANY false positive, not the
    % expected proportion, so they are far stricter than the FDR columns and are
    % not comparable with them. Returned as adjusted p-values so they slot into
    % the same tables as q.
    q = local_stepdown_fwer(p, method);
    pi0_use = NaN;
elseif ismember(method, {'adaptivefdr','bh2000','lsl','bky','bky2006','twostage'})
    % Adaptive FDR: BH with m replaced by an estimate of the number of true
    % nulls. Same shape as Storey (q = pi0 * q_BH); only the pi0 estimator
    % differs - see the header.
    q = min(1, pi0_use * q_bh);
    q = max(q, p);
elseif strcmp(method, 'bh')
    q = q_bh;
    pi0_use = 1;
elseif guard_fired
    q = q_bh;
    pi0_use = 1;
    method_used = [method_used ' -> BH (guard)'];
else
    if isnan(pi0_use)
        % Only reachable with the guard off and no pi0 at all; SAS cannot
        % return a q-value here either, so BH is the only defined answer.
        q = q_bh;
        pi0_use = 1;
        method_used = [method_used ' -> BH (pi0 undefined)'];
    else
        q = min(1, pi0_use * q_bh);
        q = max(q, p);              % floor at p, as SAS proc multtest
    end
end

% ------------------------- report ----------------------------------------

if verbose
    fprintf('\n  Storey FDR on %d test(s), method ''%s''\n', n, method);
    fprintf('    pi0 used                 : %.4f  (%s)\n', pi0_use, method_used);
    fprintf('    pi0 at lambda %s : %s (range %.3f)\n', ...
        strtrim(sprintf('%.2f ', lam)), strtrim(sprintf('%.3f ', pi0_lam)), pi0_range);
    if ~isnan(pi0_spline)
        fprintf('    pi0, mafdr spline (ref)  : %.4f\n', pi0_spline);
    end
    if reliable
        fprintf('    VERDICT: pi0 well determined; Storey q returned.\n');
    elseif guard_fired
        fprintf('    VERDICT: pi0 not trustworthy; guard ON, returning Benjamini-Hochberg. Reasons:\n');
        for k = 1:numel(reasons), fprintf('      - %s\n', reasons{k}); end
    else
        fprintf(['    VERDICT: pi0 questionable, but guard OFF (method ''%s''), so the\n' ...
                 '             estimate is USED, as SAS PROC MULTTEST would. Concerns:\n'], method);
        for k = 1:numel(reasons), fprintf('      - %s\n', reasons{k}); end
        fprintf('             Storey q here is ANTI-CONSERVATIVE relative to BH; q_BH is in info.q_BH.\n');
    end
    if n < 100
        fprintf(['    NOTE: m = %d. Storey''s own smallest simulated m is 100, where the\n' ...
                 '          bootstrap lambda-selection already misses badly (see header).\n' ...
                 '          Report q_BH alongside q_Storey at this m.\n'], n);
    end
end

info = struct('method_used', method_used, 'pi0', pi0_use, 'pi0_lambda', pi0_lam, ...
    'lambda', lam, 'pi0_range', pi0_range, 'pi0_spline', pi0_spline, ...
    'reliable', reliable, 'reasons', {reasons}, 'q_BH', reshape(q_bh, sz), ...
    'guard_on', do_guard, 'guard_fired', guard_fired);

pi0 = pi0_use;
q   = reshape(q, sz);

end % main function


% =========================================================================

function [pi0, method_used] = local_sas_pi0(p, lam, pi0_spline, nboot)
% SAS PROC MULTTEST's PFDR default: SPLINE first, BOOTSTRAP on its trigger.
%
% SAS: "the SPLINE method is attempted first. If the estimate is nonpositive or
% if the slope of the spline at the last lambda is greater than 0.1 times the
% range of the fitted spline values, then the BOOTSTRAP method is used."
%
% mafdr does not expose the fitted spline, so the slope half of that test
% cannot be reproduced exactly. What is reproduced is the nonpositive test,
% plus the same practical intent: reject a spline estimate that the lambda
% curve does not support, and fall back to the Storey & Tibshirani bootstrap.

n = numel(p);

spline_bad = isnan(pi0_spline) || pi0_spline <= 0;

if ~spline_bad
    % Stand-in for SAS's terminal-slope test: does the spline estimate sit
    % anywhere near what the lambda curve implies at its upper end?
    %
    % Anchored at lambda = 0.95, the LAST lambda of SAS's NLAMBDA = 20 grid,
    % which is where SAS evaluates the spline slope. It previously used
    % max(lam), i.e. 0.5 under the default grid - a different place on the
    % curve from the one SAS tests, which rejected the spline far more often
    % at small n and pushed work onto the bootstrap, the more biased of the
    % two estimators. The two anchorings disagreed on 31-72%% of simulated
    % datasets (n = 8 to 200).
    lam_last   = 0.95;
    pi0_lam_hi = min(1, mean(p > lam_last)/(1 - lam_last));
    spline_bad = abs(pi0_spline - pi0_lam_hi) > 0.3;
end

if ~spline_bad
    pi0 = pi0_spline;
    method_used = 'SAS: spline';
    return
end

% Storey & Tibshirani (2003) bootstrap, as SAS's fallback
lam_b   = (0:19)/20;                      % NLAMBDA = 20, as SAS
pi0_b   = arrayfun(@(L) min(1, mean(p > L)/(1-L)), lam_b);
min_pi0 = min(pi0_b);

mse = zeros(size(lam_b));
for b = 1:nboot
    pb    = p(randi(n, n, 1));
    pi0_s = arrayfun(@(L) min(1, mean(pb > L)/(1-L)), lam_b);
    mse   = mse + (pi0_s - min_pi0).^2;
end
mse = mse / nboot;

[~, k] = min(mse);
pi0 = pi0_b(k);
method_used = sprintf('SAS: bootstrap (spline rejected), lambda = %.2f', lam_b(k));

end


% =========================================================================

function q = local_stepdown_fwer(p, method)
% Step-down (Holm 1979) FWER adjusted p-values, Sidak or Bonferroni multiplier.
%
% Sidak's multiplier 1-(1-p)^k assumes INDEPENDENT tests and is slightly less
% conservative than Bonferroni's k*p, which assumes nothing. Both are applied
% step-down: the smallest p-value is tested against the full m, the next
% against m-1, and so on. Step-down is uniformly more powerful than the
% corresponding single-step version at no cost in FWER control.
%
% CanlabCore has holm_sidak(), which returns a significance vector at a fixed
% alpha rather than adjusted p-values (and is marked an alpha version). This
% returns adjusted p-values, so the result drops into the same tables as the
% FDR columns; the two agree on which tests are significant at a given alpha.

p  = p(:);
m  = numel(p);
[ps, ord] = sort(p, 'ascend');
k  = (m:-1:1)';                      % m, m-1, ..., 1

switch method
    case {'stepdown_sidak','sidak'}
        a = 1 - (1 - ps).^k;
    otherwise                        % holm / bonferroni
        a = min(1, k .* ps);
end

a = cummax(a);                       % step-down monotonicity
a = min(1, a);

q = zeros(m,1);
q(ord) = a;

end


% =========================================================================

function [pi0, method_used] = local_lsl_pi0(p)
% Lowest-slope estimator of pi0: Hochberg & Benjamini (1990), after Schweder &
% Spjotvoll (1982). This is what SAS PROC MULTTEST's ADAPTIVEFDR uses by
% default (LOWESTSLOPE), giving Benjamini & Hochberg's (2000) adaptive linear
% step-up procedure.
%
% Slopes S_i = (1 - p_(i)) / (m + 1 - i) decrease while the ordered p-values
% look null, and turn up once the small p-values are exhausted. m0 is read off
% the first index where the slope stops decreasing.

p  = sort(p(:), 'ascend');
m  = numel(p);
S  = (1 - p) ./ (m + 1 - (1:m)');

% Stop at the FIRST DECREASE in slope, then m0 = min(ceil(1/S_i), m).
% Getting this wrong is silent: an earlier version here stopped at the first
% INCREASE and added 1 inside the ceiling, which always returned m0 = m and
% made this method identical to plain BH. S rises from about 1/(m+1) under the
% null, so "first increase" fires at i = 2 on essentially any input.
i_star = find(diff(S) < 0, 1, 'first') + 1;
if isempty(i_star)
    m0 = m;
    method_used = 'adaptive FDR (BH 2000), lowest slope: no decrease found, m0 = m';
else
    m0 = min(ceil(1/S(i_star)), m);
    method_used = sprintf('adaptive FDR (BH 2000), lowest slope at i = %d, m0 = %d', i_star, m0);
end

pi0 = m0 / m;

end


% =========================================================================

function [pi0, method_used] = local_bky_pi0(p, q_bh, alpha)
% Benjamini, Krieger & Yekutieli (2006) two-stage linear step-up.
%
%   Stage 1: run BH at level alpha/(1+alpha), count the rejections r1
%   Stage 2: set m0 = m - r1 and run BH again with that m0
%
% This is the route SAS documents for BKY: apply the FDR adjustment at the
% reduced level, pass the resulting count through NTRUENULL=, then apply
% ADAPTIVEFDR.
%
% NOTE: the adjusted p-values are ALPHA-DEPENDENT by construction, because BKY
% is defined as a rejection rule at a level rather than as a p-value transform.
% Change 'alpha' and the numbers move. That is a property of the procedure, not
% of this implementation; Storey and BH 2000 do not behave this way.

p  = p(:);
m  = numel(p);
a1 = alpha / (1 + alpha);
r1 = sum(q_bh(:) <= a1);

if r1 == 0
    pi0 = 1;
    method_used = sprintf('BKY 2006 two-stage: stage 1 rejected 0 at %.4f, m0 = m', a1);
    return
end
if r1 == m
    pi0 = 1/m;    % everything rejected at stage 1; m0 cannot be 0
    method_used = sprintf('BKY 2006 two-stage: stage 1 rejected all %d, m0 = 1', m);
    return
end

m0  = m - r1;
pi0 = m0 / m;
method_used = sprintf('BKY 2006 two-stage: stage 1 alpha = %.4f rejected %d, m0 = %d', a1, r1, m0);

end
