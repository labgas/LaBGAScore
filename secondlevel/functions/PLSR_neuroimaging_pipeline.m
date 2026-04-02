function results = PLSR_neuroimaging_pipeline(X,Y,opts)

% Robust PLS regression pipeline for neuroimaging feature matrices.
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
%                             Used to compute results.Q2_global.
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
%     results.featureStability  [p x 1] proportion of runs where |beta|>0
%                        (in PLSR this is typically ~1 because betas are
%                         rarely exactly zero; included for plot symmetry)
%
%   Global baseline (interpretation only):
%     results.Q2_global   scalar  predictive-style Q2 of linear model on a
%                        global summary feature, computed as:
%                        1 - SSE / SST relative to the sample mean of Y
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
%     results.selectionFrequency [p x 1] frequency of appearing in topK
%                        absolute weights across runs (topK fixed at 20)
%
%   Permutation test:
%     results.allpermQ2     [nPerm x 1] permuted Q2 distribution
%     results.permQ2        scalar mean permuted Q2
%     results.permutation_p scalar p = mean(permQ2 >= observed Q2)
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

rng(1)

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
cvObserved  = [];
cvPredicted = [];
cvRepeatID  = [];
cvSubjectID = [];

for r = 1:opts.nRepeats

    cvOuter = cvpartition(n,'KFold',opts.outerK);

    for k = 1:opts.outerK

        trainIdx = training(cvOuter,k);
        testIdx  = test(cvOuter,k);

        ytrain = Y(trainIdx);
        ytest  = Y(testIdx);

        Xtrain = X(trainIdx,:);
        Xtest  = X(testIdx,:);

        %% leakage-free scaling (generic)
        [Xtrain, Xtest] = applyScaling(Xtrain, Xtest, opts.scale);

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

                ytr = ytrain(tr);
                yva = ytrain(va);

                if numel(ytr) < 3 || numel(yva) < 2
                    continue
                end

                [~,~,~,~,beta] = plsregress(Xtrain(tr,:),ytr,lv);

                yhat = [ones(sum(va),1) Xtrain(va,:)] * beta;

                denom = sum((yva - mean(ytr)).^2);
                if denom <= 0
                    continue
                end

                foldQ2(f) = 1 - sum((yva - yhat).^2) / denom;

            end

            innerQ2(lv) = nanmean(foldQ2);

        end

        [~,bestLV] = max(innerQ2);
        selectedLV(r,k) = bestLV;

        %% fit model

        [~,~,~,~,beta] = plsregress(Xtrain,ytrain,bestLV);
        betaStore(:,k,r) = beta;
        featureWeights(:,(r-1)*opts.outerK+k) = beta(2:end);

        yhat = [ones(sum(testIdx),1) Xtest] * beta;

        % New: collect held-out observed and predicted values
        cvObserved  = [cvObserved;  ytest(:)]; %#ok<AGROW>
        cvPredicted = [cvPredicted; yhat(:)];  %#ok<AGROW>
        cvRepeatID  = [cvRepeatID; repmat(r, sum(testIdx), 1)]; %#ok<AGROW>
        cvSubjectID = [cvSubjectID; find(testIdx)];

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
results.Q2      = nanmean(Q2(:));
results.allMSE  = MSE;
results.MSE     = nanmean(MSE(:));
results.allRMSE = RMSE;
results.RMSE    = nanmean(RMSE(:));
results.allMAE  = MAE;
results.MAE     = nanmean(MAE(:));
results.allCorr = Corr;
results.Corr    = nanmean(Corr(:));

results.selectedLV      = selectedLV;
results.betaStore       = betaStore;
results.featureWeights  = featureWeights;

results.featureStability  = mean(abs(featureWeights) > 0,2,'omitnan');
results.meanFeatureWeight = mean(featureWeights,2,'omitnan');

% New: export held-out observed/predicted values
results.cvObserved  = cvObserved;
results.cvPredicted = cvPredicted;
results.cvRepeatID = cvRepeatID;
results.cvSubjectID = cvSubjectID;

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

if std(Y) > 0 && std(yhatGlobal) > 0
    Cg = corrcoef(Y,yhatGlobal);
    results.Corr_global = Cg(1,2);
else
    results.Corr_global = NaN;
end

%% -------------------------------------------------
% 4. Final model (interpretation only)
%% -------------------------------------------------

[Xz,~] = applyScaling(X, X, opts.scale);

finalLV = round(nanmedian(selectedLV(:)));
finalLV = max(1, min(finalLV, capLV(opts.maxLV, Xz)));

[XL,YL,XS,YS,beta,PCTVAR,MSEfinal,stats] = plsregress(Xz,Y,finalLV);

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
meanBeta = nanmean(betaMat,2);
sdBeta = nanstd(betaMat,[],2);

results.meanBeta = meanBeta;
results.sdBeta = sdBeta;
results.stabilityZ = meanBeta ./ sdBeta;

%% -------------------------------------------------
% 7. Permutation testing
%% -------------------------------------------------

permQ2 = nan(opts.nPerm,1);

parfor i = 1:opts.nPerm
    yp = Y(randperm(n));
    permQ2(i) = quickCV(X,yp,opts);
end

results.allpermQ2 = permQ2;
results.permQ2 = nanmean(permQ2);
results.permutation_p = mean(permQ2 >= results.Q2,'omitnan');

figure
histogram(permQ2(~isnan(permQ2)))
hold on
xline(results.Q2)
title('Permutation Q2 distribution')

%% -------------------------------------------------
% 8. Bootstrap Q2 CI (out-of-bag bootstrap)
%% -------------------------------------------------

bootQ2 = nan(opts.nBoot,1);

parfor b = 1:opts.nBoot
    bootQ2(b) = bootstrapOOB_PLSR(X, Y, opts);
end

results.allbootQ2 = bootQ2;
results.bootQ2 = nanmean(bootQ2);
results.Q2_CI = prctile(bootQ2(~isnan(bootQ2)), [2.5 97.5]);

figure
histogram(bootQ2(~isnan(bootQ2)))
hold on
xline(results.Q2)
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

    m = sizes(i);
    idx = randsample(n,m);

    lcQ2(i) = quickCV(X(idx,:),Y(idx),opts);

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
signStability = mean(sign(results.meanFeatureWeight) == sign(betaFlat),2,'omitnan');
results.signStability = signStability;

%% ---------------------------
% 11. Top-K selection Frequency
%% ---------------------------

if size(X,2) < 20
    topK = size(X,2);
else
    topK = 20;
end

freq = zeros(size(results.meanFeatureWeight));

nRuns = size(results.featureWeights,2);

for i = 1:nRuns
    wi = results.featureWeights(:,i);
    if all(isnan(wi))
        continue
    end
    [~,idx] = sort(abs(wi),'descend');
    nKeep = min(topK, numel(idx));
    freq(idx(1:nKeep)) = freq(idx(1:nKeep)) + 1;
end

freq = freq / nRuns;
results.selectionFrequency = freq;

end

%% -------------------------------------------------
% inline / local helper functions
%% -------------------------------------------------

function maxLV = capLV(maxLVopt, Xtrain)
% capLV: ensures LV count is valid for any p/n ratio
nTr = size(Xtrain,1);

% rank-based cap is important even when p<n
rX = rank(Xtrain);

% plsregress needs lv <= min(rank(X), n-1) typically; keep your -1/-2 safeguards
maxLV = min([maxLVopt, rX-1, nTr-2]);

if isnan(maxLV) || isinf(maxLV)
    maxLV = 0;
end

maxLV = floor(maxLV);
end

function [XtrS, XteS] = applyScaling(Xtr, Xte, mode)
% leakage-free scaling options
switch lower(mode)
    case 'none'
        XtrS = Xtr;
        XteS = Xte;
    case 'center'
        mu = mean(Xtr,1);
        XtrS = Xtr - mu;
        XteS = Xte - mu;
    otherwise % 'zscore'
        mu = mean(Xtr,1);
        sd = std(Xtr,0,1);
        sd(sd==0) = 1;
        XtrS = (Xtr - mu) ./ sd;
        XteS = (Xte - mu) ./ sd;
end
end

function Q2 = quickCV(X,Y,opts)

n = length(Y);
K = min(opts.outerK,floor(n/2));

if K < 2
    Q2 = NaN;
    return
end

cv = cvpartition(n,'KFold',K);

q2 = nan(K,1);

for k = 1:K

    tr = training(cv,k);
    te = test(cv,k);

    ytr = Y(tr);
    yte = Y(te);

    if numel(ytr) < 3 || numel(yte) < 2
        continue
    end

    Xtr = X(tr,:);
    Xte = X(te,:);

    [Xtr, Xte] = applyScaling(Xtr, Xte, opts.scale);

    lv = capLV(opts.maxLV, Xtr);
    if lv < 1
        continue
    end

    [~,~,~,~,beta] = plsregress(Xtr,ytr,lv);

    yhat = [ones(sum(te),1) Xte] * beta;

    denom = sum((yte - mean(ytr)).^2);
    if denom <= 0
        continue
    end

    q2(k) = 1 - sum((yte - yhat).^2) / denom;

end

Q2 = nanmean(q2);

end

function Q2 = bootstrapOOB_PLSR(X, Y, opts)
% bootstrapOOB_PLSR
% Out-of-bag bootstrap Q2 for PLSR:
% - bootstrap sample used for training
% - OOB subjects used for testing
% - LV selected by inner CV within the bootstrap sample

n = length(Y);

% Bootstrap sample
idxBoot = randsample(n, n, true);

% OOB = subjects not selected at least once
inBag = false(n,1);
inBag(idxBoot) = true;
oob = ~inBag;

% Need enough OOB cases
if sum(oob) < 2
    Q2 = NaN;
    return
end

ytrain = Y(idxBoot);
ytest  = Y(oob);

if numel(ytrain) < 4 || numel(ytest) < 2
    Q2 = NaN;
    return
end

Xtrain = X(idxBoot,:);
Xtest  = X(oob,:);

% leakage-free scaling
[Xtrain, Xtest] = applyScaling(Xtrain, Xtest, opts.scale);

% cap LV
maxLV = capLV(opts.maxLV, Xtrain);
if maxLV < 1
    Q2 = NaN;
    return
end

% inner CV for LV tuning
innerK = min([opts.innerK, length(ytrain)-1]);
if innerK < 2
    Q2 = NaN;
    return
end

cvInner = cvpartition(length(ytrain),'KFold',innerK);

innerQ2 = nan(maxLV,1);

for lv = 1:maxLV
    foldQ2 = nan(innerK,1);

    for f = 1:innerK
        tr = training(cvInner,f);
        va = test(cvInner,f);

        ytr = ytrain(tr);
        yva = ytrain(va);

        if numel(ytr) < 3 || numel(yva) < 2
            continue
        end

        [~,~,~,~,beta] = plsregress(Xtrain(tr,:), ytr, lv);
        yhat = [ones(sum(va),1) Xtrain(va,:)] * beta;

        denom = sum((yva - mean(ytr)).^2);
        if denom <= 0
            continue
        end

        foldQ2(f) = 1 - sum((yva - yhat).^2) / denom;
    end

    innerQ2(lv) = nanmean(foldQ2);
end

[~,bestLV] = max(innerQ2);

% final fit on bootstrap sample
[~,~,~,~,beta] = plsregress(Xtrain, ytrain, bestLV);

yhat = [ones(sum(oob),1) Xtest] * beta;
denom = sum((ytest - mean(ytrain)).^2);

if denom <= 0
    Q2 = NaN;
else
    Q2 = 1 - sum((ytest - yhat).^2) / denom;
end

end