%% control_mvpa_reg_cov_permutation.m
%
% POSITIVE CONTROL for the permutation test added to prep_3a's
% domvpa_reg_cov block (prep_3a v_next, 2026/09/23).
%
% WHY THIS EXISTS
% predict() returns pred_outcome_r, mse, rmse, meanabserr and cverr, and
% nothing inferential - there is no p-value anywhere in its output struct.
% The permutation block supplies one. A test of that block has to show two
% things, and the second is the one that actually matters:
%
%   1. POWER   - with real signal, the test detects it.
%   2. CALIBRATION - with NO signal, p is uniform and the null is centred on
%                    zero. A test that only passes (1) can still be useless:
%                    an anti-conservative null would "detect" signal that is
%                    not there, which is the failure mode worth guarding.
%
% Fully synthetic, so the truth is known exactly rather than assumed:
% X is Gaussian noise; under the null Y is independent Gaussian noise; under
% the alternative Y is a sparse linear function of X plus noise.
%
% EXPECTED
%   NULL      mean pred_outcome_r ~ 0, permutation null centred on ~0,
%             p-values approximately uniform, ~5% below .05
%   SIGNAL    pred_outcome_r clearly positive, p at the permutation floor
%   BOTH      the permutation null is centred near zero. It is built by
%             permuting Y, so it must not shift with real signal - if it
%             does, the p-values are wrong.
%
% If a change to the permutation block makes this control fail, the block is
% broken.
%
% -------------------------------------------------------------------------
%
% part of the LaBGAScore script-checker positive controls

rng(20260923, 'twister');

nsub    = 60;      % subjects
nvox    = 500;     % voxels
nfolds  = 5;
nperm   = 200;     % permutations per dataset
ndataset_null = 20;

fprintf('\n=== MVPA covariate permutation test: positive control ===\n');
fprintf('n = %d subjects, %d voxels, %d-fold CV, %d permutations\n\n', ...
    nsub, nvox, nfolds, nperm);

% ---- 1. CALIBRATION: no signal anywhere ---------------------------------
fprintf('--- NULL datasets (Y independent of X) ---\n');
pn = nan(ndataset_null,1); rn = nan(ndataset_null,1); nullmean = nan(ndataset_null,1);
for ds = 1:ndataset_null
    X = randn(nvox, nsub);
    Y = randn(nsub, 1);
    [rn(ds), pn(ds), nv] = run_one(X, Y, nfolds, nperm);
    nullmean(ds) = mean(nv);
    fprintf('  dataset %2d: r = %+.4f, p = %.4f, null mean = %+.4f\n', ds, rn(ds), pn(ds), nullmean(ds));
end

fprintf('\n  observed r  : mean %+.4f (should be ~0, slightly negative is normal for CV)\n', mean(rn));
fprintf('  null centre : mean of null means %+.4f (MUST be ~0)\n', mean(nullmean));
fprintf('  p-values    : mean %.3f (uniform => ~0.5), %d of %d below .05 (expect ~%.1f)\n', ...
    mean(pn), sum(pn < .05), ndataset_null, 0.05*ndataset_null);

ok_centre  = abs(mean(nullmean)) < 0.05;
ok_uniform = sum(pn < .05) <= max(3, ceil(0.05*ndataset_null*3));

% ---- 2. POWER: real signal ----------------------------------------------
fprintf('\n--- SIGNAL dataset (Y is a sparse linear function of X) ---\n');
% The signal has to live in a HIGH-VARIANCE direction of X, because that is
% what principal-component regression can actually use. An earlier version of
% this control put the signal in 50 of 500 isotropic Gaussian voxels: the
% truth was exactly known, but it was spread evenly across all components, so
% cv_pcr could not recover it and returned r = -0.24. That tested the
% synthetic design, not the permutation machinery.
%
% Real fMRI contrast data has strong spatial covariance and the effects of
% interest usually do project onto leading components. So: one latent score
% z drives both a spatial pattern in X and the outcome Y.
z       = randn(nsub, 1);                       % latent subject score
pattern = randn(nvox, 1);                       % its spatial expression
X = pattern * (z' * 3) + randn(nvox, nsub);     % voxels x subjects, pattern dominates
Y = z + randn(nsub,1) * 1.0;                    % R^2 = 0.5 against the latent score
[rs, ps, nvs] = run_one(X, Y, nfolds, nperm);
fprintf('  r = %+.4f, p = %.4f, null mean = %+.4f\n', rs, ps, mean(nvs));

ok_power   = rs > 0.3 && ps < 0.01;
ok_centre2 = abs(mean(nvs)) < 0.05;

% ---- verdict ------------------------------------------------------------
fprintf('\n=== VERDICT ===\n');
fprintf('  null centred on zero (no-signal)   : %s\n', pass(ok_centre));
fprintf('  p approximately uniform            : %s\n', pass(ok_uniform));
fprintf('  detects real signal                : %s\n', pass(ok_power));
fprintf('  null STILL centred with signal     : %s\n', pass(ok_centre2));
if ok_centre && ok_uniform && ok_power && ok_centre2
    fprintf('\nPERMUTATION TEST OK - calibrated under the null and powered under signal.\n');
else
    fprintf('\nPERMUTATION TEST FAILED - see the lines above.\n');
end

function s = pass(tf)
    if tf, s = 'PASS'; else, s = '*** FAIL ***'; end
end


function [obs_r, p_perm, nullv] = run_one(X, Y, nfolds, nperm)
% Reproduces prep_3a's domvpa_reg_cov permutation logic exactly:
% seeded folds, predict() on the real outcome, then the same folds
% re-used while the outcome is permuted and the whole CV re-run.
[nvox, nsub] = size(X);
    d = fmri_data;
    d.dat = X;                 % voxels x subjects
    d.Y   = Y(:);
    d.removed_voxels = false(nvox,1);
    d.removed_images = false(nsub,1);

    cv = cvpartition(nsub, 'KFold', nfolds);
    fold_labels = zeros(nsub,1);
    for k = 1:cv.NumTestSets, fold_labels(cv.test(k)) = k; end

    [~, st] = predict(d, 'algorithm_name', 'cv_pcr', 'nfolds', fold_labels, ...
                      'error_type', 'mse', 'verbose', 0);
    obs_r = corr(st.yfit, d.Y);

    % permutations: permute Y, re-run the WHOLE CV, folds held fixed
    permidx = zeros(nsub, nperm);
    for pp = 1:nperm, permidx(:,pp) = randperm(nsub)'; end
    nullv = nan(nperm,1);
    for pp = 1:nperm
        dp = d; dp.Y = d.Y(permidx(:,pp));
        [~, sp] = predict(dp, 'algorithm_name', 'cv_pcr', 'nfolds', fold_labels, ...
                          'error_type', 'mse', 'verbose', 0);
        nullv(pp) = corr(sp.yfit, dp.Y);
    end
    p_perm = (sum(nullv >= obs_r) + 1) / (nperm + 1);
end
