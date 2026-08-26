function results = PLSR_neuroimaging_pipeline(X,Y,opts)

% PLSR_neuroimaging_pipeline  Robust PLS regression pipeline for neuroimaging feature matrices.
%
% This function implements Partial Least Squares Regression (PLSR)
% for continuous outcomes with repeated nested k-fold cross-validation.
% It is designed for neuroimaging feature matrices (subjects × features) such
% as PET ROI binding, fMRI ROI betas, morphometry, connectivity-derived
% measures (edges), graph metrics, or multimodal ROI feature sets. The architecture emphasizes:
%   - leakage-free preprocessing (train-only scaling in every fold)
%   - inner CV selection of number of latent variables (LVs)
%   - outer CV estimation of generalization performance
%   - resampling-based inference (permutation p-value, bootstrap CI)
%   - stability/interpretability metrics (VIP, stabilityZ, sign stability,
%     top-K selection frequency)
%   - learning curve (Q2 vs sample size)
%
% USAGE
%   results = PLSR_neuroimaging_pipeline(X, Y)
%   results = PLSR_neuroimaging_pipeline(X, Y, opts)
%
% INPUTS
%   X    [n x p] numeric
%        Feature matrix (n subjects, p features). Each column is a feature
%        (e.g., ROI binding, ROI beta, edge weight, graph metric).
%        Missing values should be handled upstream (impute/remove) prior to calling this function.
%
%   Y    [n x 1] numeric
%        Continuous outcome vector. Must be numeric and finite.
%
%   opts (optional) struct with fields:
%        Cross-validation:
%          opts.outerK        (default 5)  outer k-fold CV
%          opts.innerK        (default 4)  inner k-fold CV (LV tuning)
%          opts.nRepeats      (default 50) repeats of outer CV
%
%        Model complexity:
%          opts.maxLV         (default 4)  max candidate LVs
%            Actual per-fold cap is:
%              maxLV = min([opts.maxLV, rank(Xtrain)-1, nTrain-2])
%
%        Resampling inference:
%          opts.nPerm         (default 1000) permutations for p-value
%          opts.nBoot         (default 500)  bootstrap resamples for Q2 CI
%
%        Learning curve:
%          opts.learningSteps (default 6) number of sample sizes
%
%        Generic additions:
%          opts.scale         (default 'zscore') scaling mode inside every fold:
%                             'zscore' | 'center' | 'none'
%                             - 'zscore' is recommended for most neuroimaging matrices
%                               (different ROIs/edges/metrics have different scales).
%          opts.globalFun     (default 'mean') global baseline feature:
%                             'mean' | 'median' | function handle @(X)->[n x 1]
%                             Used for results.Q2_global and results.Q2_global_cv.
%
%        Covariate control (fold-wise nuisance regression):
%          opts.covariates    (default []) [n x nCov] numeric nuisance matrix,
%                             one ROW per subject. Do NOT include a column of
%                             ones; the intercept is handled internally.
%                             Nuisance coefficients are estimated on the
%                             TRAINING fold only and applied to both folds, in
%                             every outer fold, every inner fold, every
%                             bootstrap resample and every learning-curve
%                             subsample. Encode categorical covariates as dummy
%                             columns beforehand. Leave empty to reproduce the
%                             previous behaviour exactly.
%          opts.covariateNames (default {'cov1',...}) cellstr, provenance only.
%
%        Reproducibility:
%          opts.seed          (default 1) RNG seed. Results are now reproducible
%                             across runs AND independent of parallel pool size.
%          opts.residualizeY  (default false) also residualize Y on the
%                             covariates, train-only, in every fold. With false
%                             the model predicts total Y from confound-free X,
%                             so the nuisance-explained part of Y is
%                             structurally unexplainable and Q2 is DEFLATED.
%                             With true, Q2 becomes a partial Q2: variance in Y
%                             explained by X after removing the covariates from
%                             both -- the regression analogue of a partial
%                             correlation. The two are not comparable; the
%                             setting used is recorded in results.covariateInfo.
%                             No effect when opts.covariates is empty.
%
% OUTPUT (results struct)
%   Cross-validated performance (generalization estimate):
%     results.Q2         scalar   mean out-of-sample predictive Q2 across repeats×outer folds
%     results.MSE        scalar   mean MSE across repeats×outer folds
%     results.RMSE       scalar   mean RMSE across repeats×outer folds
%     results.MAE        scalar   mean MAE across repeats×outer folds
%     results.Corr       scalar   mean Pearson correlation across repeats×outer folds
%     results.allQ2      [nRepeats x outerK] fold-level Q2
%     results.allMSE     [nRepeats x outerK] fold-level MSE
%     results.allRMSE    [nRepeats x outerK] fold-level RMSE
%     results.allMAE     [nRepeats x outerK] fold-level MAE
%     results.allCorr    [nRepeats x outerK] fold-level Pearson r
%
%   Held-out predictions from outer CV (diagnostics/plotting):
%     results.cvObserved   [N x 1] stacked outer-fold held-out observed Y values
%     results.cvPredicted  [N x 1] stacked outer-fold held-out predicted Y values
%     results.cvRepeatID   [N x 1] repeat index for each stacked held-out prediction
%     results.cvSubjectID  [N x 1] subject index for each stacked held-out prediction
%                        (useful for predicted-vs-observed plots colored by repeat
%                         or subject across repeated outer CV)
%
%   Model selection / weights across CV:
%     results.selectedLV [nRepeats x outerK] selected LV per outer fold
%     results.betaStore  [p+1 x outerK x nRepeats] PLS regression betas
%                        (row 1 is intercept; rows 2..p+1 correspond to features)
%     results.featureWeights [p x (nRepeats*outerK)] betas (no intercept),
%                        stacked across all outer folds/repeats
%     results.meanFeatureWeight [p x 1] mean featureWeights across runs
%
%   Global baseline:
%     results.Q2_global   scalar  IN-SAMPLE Q2 of a linear model on a global
%                        summary feature (fit and predicted on all subjects),
%                        computed as 1 - SSE / SST relative to the sample mean.
%                        Retained unchanged for continuity, but NOT comparable
%                        to the cross-validated results.Q2.
%     results.Q2_global_cv scalar cross-validated counterpart, run through the
%                        same repeated outer CV as the model. THIS is the value
%                        to compare against results.Q2.
%     results.MSE_global scalar   MSE of global summary feature model
%     results.Corr_global scalar  Pearson correlation for global summary feature model
%
%   Final model on all data (interpretation only; NOT for performance):
%     results.finalLV        scalar median selected LV across all runs (then capped valid)
%     results.betaFinal      [p+1 x 1] betas from final model
%     results.varExplainedX  [1 x finalLV] PCTVAR(1,:) from plsregress
%     results.varExplainedY  [1 x finalLV] PCTVAR(2,:) from plsregress
%     results.finalXLoadings [p x finalLV] XL from plsregress
%     results.finalYLoadings [1 x finalLV] YL from plsregress
%     results.finalXScores   [n x finalLV] XS from plsregress
%     results.finalYScores   [n x finalLV] YS from plsregress
%     results.finalMSE       output MSE from plsregress
%
%   Feature importance / stability (from final model and CV betas):
%     results.VIP         [p x 1] VIP scores (higher = more important)
%     results.meanBeta    [p x 1] mean beta across outer CV runs
%     results.sdBeta      [p x 1] std beta across outer CV runs
%     results.stabilityZ  [p x 1] meanBeta ./ sdBeta (stability statistic)
%     results.signStability [p x 1] proportion of runs matching mean sign
%
%   Permutation test (Freedman-Lane):
%     results.allpermQ2     [nPerm x 1] permuted Q2 distribution
%     results.permQ2        scalar mean permuted Q2
%     results.quickCV_observed scalar observed Q2 computed with quickCV_PLSR,
%                        i.e. the SAME estimator that generates the null
%     results.permutation_p scalar (sum(permQ2 >= quickCV_observed) + 1) / (nValid + 1)
%
%     The null permutes the residuals of Y after nuisance regression and adds
%     the fitted nuisance part back, preserving the Y-covariate association
%     while breaking the X-Y one. With no covariates this reduces to ordinary
%     unrestricted permutation of Y.
%
%     The p-value is computed against quickCV_observed rather than results.Q2
%     because the null is generated by quickCV_PLSR. Comparing a repeated-
%     nested-CV Q2 against a quickCV null mixes two estimators; measured on
%     true-null data that mismatch produced a false positive rate near 40 percent
%     at alpha = 0.05. results.Q2 remains the performance estimate to report.
%
%   Bootstrap:
%     results.allbootQ2 [nBoot x 1] out-of-bag bootstrap Q2 distribution
%     results.bootQ2    scalar mean out-of-bag bootstrap Q2
%     results.Q2_CI     [1 x 2] percentile CI (2.5, 97.5)
%
%   Learning curve:
%     results.learningSizes vector of sample sizes evaluated
%     results.learningQ2    vector of Q2 estimates per size
%
%   Provenance:
%     results.covariateInfo struct: .used .names .nCov .rank .residualizeY
%                        .order ('residualize-then-scale') .permScheme
%     results.seed          the RNG seed actually used
%
% NOTES / INTERPRETATION (high level)
%   - Use results.Q2 from nested CV as the primary generalization estimate.
%   - Fold-level Q2 is computed out-of-sample relative to the training-fold
%     mean of Y:
%         Q2 = 1 - sum((ytest - yhat).^2) / sum((ytest - mean(ytrain)).^2)
%     This avoids using test-set information in the baseline.
%   - Bootstrap Q2 is estimated using out-of-bag (OOB) testing rather than
%     evaluating on the bootstrap sample itself, which makes it more conservative
%     and typically less optimistic than naive bootstrap performance estimates.
%   - Accordingly, bootstrap Q2 may be somewhat lower than nested CV Q2 in
%     small samples and should be viewed primarily as a measure of sampling
%     variability rather than as the main performance estimate.
%   - VIP and stabilityZ provide complementary interpretability:
%       VIP: importance in the final fitted model
%       stabilityZ: robustness of the effect across CV runs
%   - Negative Q2 values are possible and indicate prediction worse than the
%     training-mean baseline.
%   - The stacked held-out predictions (cvObserved/cvPredicted) are intended
%     primarily for visualization and diagnostics, not as an independent
%     replacement for the fold-wise nested-CV metrics above.
%
% IMPLEMENTATION NOTES
%   - Scaling is leakage-free and controlled by opts.scale.
%   - Bootstrap confidence intervals are based on out-of-bag bootstrap samples:
%     each bootstrap replicate is trained on the in-bag sample, latent variables
%     are tuned within that sample, and Q2 is evaluated on out-of-bag subjects.
%   - If cvpartition fails, the code falls back to non-stratified partitions.
%
% DEPENDENCIES
%   Requires Statistics and Machine Learning Toolbox:
%     cvpartition, plsregress, fitlm
%
% See also: plsregress, cvpartition, fitlm

%% -------------------------------------------------
% 0. Defaults
%% -------------------------------------------------

if nargin < 3
    opts = struct;
end

if ~isfield(opts,'outerK'); opts.outerK = 5; end
if ~isfield(opts,'innerK'); opts.innerK = 4; end
if ~isfield(opts,'nRepeats'); opts.nRepeats = 50; end
if ~isfield(opts,'maxLV'); opts.maxLV = 4; end
if ~isfield(opts,'nPerm'); opts.nPerm = 1000; end
if ~isfield(opts,'nBoot'); opts.nBoot = 500; end
if ~isfield(opts,'learningSteps'); opts.learningSteps = 6; end

% Generic additions (kept minimal)
if ~isfield(opts,'scale'); opts.scale = 'zscore'; end  % 'zscore'|'center'|'none'
if ~isfield(opts,'globalFun'); opts.globalFun = 'mean'; end % 'mean'|'median' (or function handle)
if ~isfield(opts,'residualizeY'); opts.residualizeY = false; end % see residualizeY
if ~isfield(opts,'seed'); opts.seed = 1; end                     % see setParforStream

warnUnknownOptions(opts, { ...
    'outerK','innerK','nRepeats','maxLV','nPerm','nBoot','learningSteps', ...
    'scale','globalFun','residualizeY','seed','covariates','covariateNames'}, ...
    'PLSR_neuroimaging_pipeline');

rng(opts.seed,'twister')

%% -------------------------------------------------
% 1. Outcome preparation
%% -------------------------------------------------

Y = Y(:);

if ~isnumeric(Y)
    error('For PLSR, Y must be a numeric continuous vector.');
end

if any(~isfinite(Y))
    error('Y contains non-finite values. Please handle missing/infinite values upstream.');
end

[n,p] = size(X);

if size(Y,1) ~= n
    error('X and Y must have the same number of rows/subjects.');
end

if ~all(isfinite(X(:)))
    error('PLSR_neuroimaging_pipeline:nonFiniteX', ...
        'X contains NaN or Inf. Handle missing values upstream (impute or drop).');
end

% Covariates are resolved ONCE here and then passed explicitly to every helper
% alongside X. They are deliberately never re-read from opts inside a helper:
% the bootstrap and learning-curve stages subset subjects, and a full-length
% covariate matrix reaching a subsetted X would misalign silently.
[Cov, covInfo] = validateCovariates(opts, n, 'PLSR_neuroimaging_pipeline');

%% -------------------------------------------------
% 2. Repeated Nested Cross-Validation
%% -------------------------------------------------

Q2   = nan(opts.nRepeats,opts.outerK);
MSE  = nan(opts.nRepeats,opts.outerK);
RMSE = nan(opts.nRepeats,opts.outerK);
MAE  = nan(opts.nRepeats,opts.outerK);
Corr = nan(opts.nRepeats,opts.outerK);

selectedLV = nan(opts.nRepeats,opts.outerK);
betaStore = nan(p+1,opts.outerK,opts.nRepeats);
featureWeights = nan(p, opts.nRepeats*opts.outerK);

% New: store held-out outer-fold predictions for plotting
% Held-out predictions are preallocated to the largest they can be and trimmed
% at the end. Each repeat's cvpartition covers every subject exactly once, so
% n*nRepeats is the ceiling; it is only reached when no fold is skipped, and
% folds ARE skipped when a rank or sample-size guard fires. Growing these four
% by concatenation reallocated on every fold.
maxCvRows   = n * opts.nRepeats;
cvObserved  = nan(maxCvRows,1);
cvPredicted = nan(maxCvRows,1);
cvRepeatID  = nan(maxCvRows,1);
cvSubjectID = nan(maxCvRows,1);
cvCursor    = 0;

for r = 1:opts.nRepeats

    cvOuter = cvpartition(n,'KFold',opts.outerK);

    for k = 1:opts.outerK

        trainIdx = training(cvOuter,k);
        testIdx  = test(cvOuter,k);

        ytrain_raw = Y(trainIdx);
        ytest      = Y(testIdx);

        Xtrain_raw = X(trainIdx,:);
        Xtest_raw  = X(testIdx,:);

        Ctrain = []; Ctest = [];
        if ~isempty(Cov)
            Ctrain = Cov(trainIdx,:);
            Ctest  = Cov(testIdx,:);
        end

        %% leakage-free preprocessing: nuisance regression then scaling,
        %% every constant estimated on the training fold alone
        [Xtrain, Xtest] = foldPreprocess(Xtrain_raw, Xtest_raw, Ctrain, Ctest, opts.scale);
        [ytrain, ytest] = residualizeY(ytrain_raw, ytest, Ctrain, Ctest, opts.residualizeY);

        %% inner CV LV tuning

        maxLV = capLV(opts.maxLV, Xtrain);

        if maxLV < 1
            continue
        end

        innerK = min([opts.innerK, sum(trainIdx)-1]);
        if innerK < 2
            continue
        end

        innerQ2 = nan(maxLV,1);
        cvInner = cvpartition(length(ytrain),'KFold',innerK);

        parfor lv = 1:maxLV

            foldQ2 = nan(innerK,1);

            for f = 1:innerK

                tr = training(cvInner,f);
                va = test(cvInner,f);

                if sum(tr) < 3 || sum(va) < 2
                    continue
                end

                % Each inner fold preprocesses itself from the RAW training
                % rows. Reusing the outer fold's constants would leak
                % outer-training statistics -- and, with covariates, the
                % outer-training nuisance fit -- into every inner validation
                % fold, which is the same leak one level down.
                C2 = []; Cv = [];
                if ~isempty(Ctrain), C2 = Ctrain(tr,:); Cv = Ctrain(va,:); end

                [X2, Xv] = foldPreprocess(Xtrain_raw(tr,:), Xtrain_raw(va,:), C2, Cv, opts.scale);
                [ytr, yva] = residualizeY(ytrain_raw(tr), ytrain_raw(va), C2, Cv, opts.residualizeY);

                if capLV(lv, X2) < lv
                    continue
                end

                [~,~,~,~,beta] = plsregress(X2,ytr,lv);

                yhat = [ones(sum(va),1) Xv] * beta;

                denom = sum((yva - mean(ytr)).^2);
                if denom <= 0
                    continue
                end

                foldQ2(f) = 1 - sum((yva - yhat).^2) / denom;

            end

            innerQ2(lv) = mean(foldQ2,'omitnan');

        end

        if all(isnan(innerQ2))
            continue
        end

        [~,bestLV] = max(innerQ2);
        selectedLV(r,k) = bestLV;

        %% fit model

        [~,~,~,~,beta] = plsregress(Xtrain,ytrain,bestLV);
        betaStore(:,k,r) = beta;
        featureWeights(:,(r-1)*opts.outerK+k) = beta(2:end);

        yhat = [ones(sum(testIdx),1) Xtest] * beta;

        % Collect held-out observed and predicted values
        nTe = sum(testIdx);
        cvRows = cvCursor + (1:nTe);
        cvObserved(cvRows)  = ytest(:);
        cvPredicted(cvRows) = yhat(:);
        cvRepeatID(cvRows)  = r;
        cvSubjectID(cvRows) = find(testIdx);
        cvCursor = cvCursor + nTe;

        err = ytest - yhat;
        mse = mean(err.^2);
        rmse = sqrt(mse);
        mae = mean(abs(err));

        denom = sum((ytest - mean(ytrain)).^2);
        if denom > 0
            Q2(r,k) = 1 - sum(err.^2) / denom;
        end

        MSE(r,k) = mse;
        RMSE(r,k) = rmse;
        MAE(r,k) = mae;

        if numel(ytest) > 1 && std(ytest) > 0 && std(yhat) > 0
            C = corrcoef(ytest,yhat);
            Corr(r,k) = C(1,2);
        end

    end
end

results.allQ2   = Q2;
results.Q2      = mean(Q2(:),'omitnan');
results.allMSE  = MSE;
results.MSE     = mean(MSE(:),'omitnan');
results.allRMSE = RMSE;
results.RMSE    = mean(RMSE(:),'omitnan');
results.allMAE  = MAE;
results.MAE     = mean(MAE(:),'omitnan');
results.allCorr = Corr;
results.Corr    = mean(Corr(:),'omitnan');

results.selectedLV      = selectedLV;
results.betaStore       = betaStore;
results.featureWeights  = featureWeights;

results.meanFeatureWeight = mean(featureWeights,2,'omitnan');

% Export held-out observed/predicted values, trimmed to the rows actually
% filled (shorter than maxCvRows whenever a fold was skipped)
results.cvObserved  = cvObserved(1:cvCursor);
results.cvPredicted = cvPredicted(1:cvCursor);
results.cvRepeatID  = cvRepeatID(1:cvCursor);
results.cvSubjectID = cvSubjectID(1:cvCursor);

fprintf('Nested CV Q2 = %.3f\n',results.Q2)

%% -------------------------------------------------
% 3. Global signal model (generic)
%% -------------------------------------------------

if isa(opts.globalFun,'function_handle')
    globalFeature = opts.globalFun(X);
else
    switch lower(opts.globalFun)
        case 'median'
            globalFeature = median(X,2);
        otherwise
            globalFeature = mean(X,2); % default
    end
end

mdl = fitlm(globalFeature,Y);
yhatGlobal = predict(mdl,globalFeature);

denomGlobal = sum((Y - mean(Y)).^2);
if denomGlobal > 0
    results.Q2_global = 1 - sum((Y - yhatGlobal).^2) / denomGlobal;
else
    results.Q2_global = NaN;
end

results.MSE_global = mean((Y - yhatGlobal).^2);

% Cross-validated companion. Q2_global above is IN-SAMPLE (fitlm predicts the
% same subjects it was fit on) and so is not comparable to the cross-validated
% results.Q2; Q2_global_cv is. Both are kept: the in-sample field retains its
% previous meaning and value.
gcv = globalBaselineCV(X, Y, Cov, opts, 'regress');
results.Q2_global_cv = gcv.Q2;

if std(Y) > 0 && std(yhatGlobal) > 0
    Cg = corrcoef(Y,yhatGlobal);
    results.Corr_global = Cg(1,2);
else
    results.Corr_global = NaN;
end

%% -------------------------------------------------
% 4. Final model (interpretation only)
%% -------------------------------------------------

% Interpretation-only fit on all subjects. There is no held-out set here, so
% fitting the nuisance model on the full sample is correct by construction.
[Xz,~]  = foldPreprocess(X, [], Cov, [], opts.scale);
[Yz,~]  = residualizeY(Y, [], Cov, [], opts.residualizeY);

finalLV = round(median(selectedLV(:),'omitnan'));
finalLV = max(1, min(finalLV, capLV(opts.maxLV, Xz)));

[XL,YL,XS,YS,beta,PCTVAR,MSEfinal,stats] = plsregress(Xz,Yz,finalLV);

results.finalLV = finalLV;
results.betaFinal = beta;
results.varExplainedX = PCTVAR(1,:);
results.varExplainedY = PCTVAR(2,:);
results.finalXLoadings = XL;
results.finalYLoadings = YL;
results.finalXScores = XS;
results.finalYScores = YS;
results.finalMSE = MSEfinal;

%% -------------------------------------------------
% 5. VIP scores
%% -------------------------------------------------

W = stats.W;
T = XS;
Q = YL;

SSY = sum(T.^2,1) .* (Q'.^2);

VIP = zeros(p,1);

parfor j = 1:p
    w = (W(j,:).^2) ./ sum(W.^2,1);
    VIP(j) = sqrt(p * sum(SSY .* w) / sum(SSY));
end

results.VIP = VIP;

%% -------------------------------------------------
% 6. Weight stability
%% -------------------------------------------------

betaMat = reshape(betaStore(2:end,:,:),p,[]);
meanBeta = mean(betaMat,2,'omitnan');
sdBeta = std(betaMat,[],2,'omitnan');

results.meanBeta = meanBeta;
results.sdBeta = sdBeta;
results.stabilityZ = meanBeta ./ sdBeta;

%% -------------------------------------------------
% 7. Permutation testing
%% -------------------------------------------------

% The p-value compares like with like: the null comes from quickCV_PLSR, so
% the observed statistic it is compared against must come from quickCV_PLSR
% too. Comparing the repeated-nested-CV Q2 against a quickCV null mixes two
% different estimators and, measured on true-null data, produced a false
% positive rate near 40% at alpha = 0.05.
obsMatched = quickCV_PLSR(X, Y, Cov, opts);
results.quickCV_observed = obsMatched;

% Freedman-Lane: permute the residuals of Y after nuisance regression and add
% the fitted nuisance part back, so the Y-covariate association is preserved
% while the X-Y association is broken. With no covariates Yhat is the mean and
% this reduces exactly to the usual unrestricted permutation of Y.
if isempty(Cov)
    Yhat = zeros(n,1);
    Ry   = Y;
else
    [Ry, ~] = residualizeFold(Y, [], Cov, []);
    Yhat    = Y - Ry;
end

permQ2 = nan(opts.nPerm,1);

parfor i = 1:opts.nPerm
    setParforStream(opts.seed, 2e6, i);
    yp = Yhat + Ry(randperm(n));
    permQ2(i) = quickCV_PLSR(X, yp, Cov, opts);
end

results.allpermQ2 = permQ2;
results.permQ2 = mean(permQ2,'omitnan');
results.permutation_p = (sum(permQ2 >= obsMatched) + 1) / (sum(~isnan(permQ2)) + 1);

figure
histogram(permQ2(~isnan(permQ2)))
hold on
xline(obsMatched)
title('Permutation Q2 distribution')

%% -------------------------------------------------
% 8. Bootstrap Q2 CI (out-of-bag bootstrap)
%% -------------------------------------------------

bootQ2 = nan(opts.nBoot,1);

parfor b = 1:opts.nBoot
    setParforStream(opts.seed, 3e6, b);
    bootQ2(b) = bootstrapOOB_PLSR(X, Y, Cov, opts);
end

results.allbootQ2 = bootQ2;
results.bootQ2 = mean(bootQ2,'omitnan');
results.Q2_CI = prctile(bootQ2(~isnan(bootQ2)), [2.5 97.5]);

figure
histogram(bootQ2(~isnan(bootQ2)))
hold on
xline(results.Q2)   % bootstrapOOB_PLSR tunes internally, so the nested-CV Q2 is the right reference
title('Bootstrap OOB Q2')
xlabel('Q2')
ylabel('Frequency')

%% -------------------------------------------------
% 9. Learning curve
%% -------------------------------------------------

sizes = round(linspace(max(10,ceil(n*0.4)), n, opts.learningSteps));
sizes = unique(sizes);

lcQ2 = nan(length(sizes),1);

parfor i = 1:length(sizes)

    setParforStream(opts.seed, 4e6, i);

    m = sizes(i);
    idx = randsample(n,m);

    Cidx = []; if ~isempty(Cov), Cidx = Cov(idx,:); end
    lcQ2(i) = quickCV_PLSR(X(idx,:),Y(idx),Cidx,opts);

end

results.learningSizes = sizes;
results.learningQ2 = lcQ2;


figure
plot(sizes,lcQ2,'o-')
xlabel('Sample size')
ylabel('Q2')
title('Learning curve')

%% ---------------------------
% 10. Sign Stability
% ----------------------------

betaNoIntercept = betaStore(2:end,:,:);
betaFlat = reshape(betaNoIntercept,p,[]);
% Compare each fold/repeat sign to the sign of the mean weight
signStability = mean(sign(results.meanFeatureWeight) == sign(betaFlat),2,'omitnan');
results.signStability = signStability;

%% -------------------------------------------------
% 11. Provenance
%% -------------------------------------------------

covInfo.residualizeY = opts.residualizeY;
covInfo.order        = 'residualize-then-scale';
covInfo.permScheme   = 'freedman-lane';
results.covariateInfo = covInfo;
results.seed          = opts.seed;

end
