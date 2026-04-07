function results = ENet_neuroimaging_pipeline(X,Y,opts)
% 
% Robust Elastic Net pipeline for neuroimaging feature matrices.
%
% This function implements Elastic Net regularized logistic-style classification
% (via lassoglm / elastic net on a binary outcome) with repeated nested k-fold
% cross-validation. It is designed for neuroimaging feature matrices
% (subjects × features) such as PET ROI binding, fMRI ROI betas, morphometry,
% connectivity-derived measures (edges), graph metrics, or multimodal ROI
% feature sets. 
% 
% The architecture emphasizes:
%   - leakage-free preprocessing (train-only scaling in every fold)
%   - inner CV selection of Elastic Net hyperparameters (alpha, lambda)
%   - outer CV estimation of generalization performance
%   - resampling-based inference (permutation p-value, bootstrap CI)
%   - stability/interpretability metrics (mean weights, non-zero stability,
%     sign stability, top-K selection frequency)
%   - learning curve (AUC vs sample size)
%
% USAGE
%   results = ENet_neuroimaging_pipeline(X, Y)
%   results = ENet_neuroimaging_pipeline(X, Y, opts)
%
% INPUTS
%   X    [n x p] numeric
%        Feature matrix (n subjects, p features). Each column is a feature
%        (e.g., ROI binding, ROI beta, edge weight, graph metric).
%        Missing values should be handled upstream (impute/remove) prior to calling this function.
%
%   Y    [n x 1] labels
%        Binary outcome labels. Can be numeric, logical, categorical, string,
%        or cellstr. Internally, Y is converted to numeric group indices (if
%        needed), then binarized as:
%            yNum = double(Y == max(Y))
%        Therefore, the "positive class" corresponds to the maximum label.
%
%   opts (optional) struct with fields:
%        Cross-validation:
%          opts.outerK        (default 5)  outer k-fold CV
%          opts.innerK        (default 4)  inner k-fold CV (hyperparameter tuning)
%          opts.nRepeats      (default 50) repeats of outer CV
%
%        Elastic Net hyperparameters:
%          opts.alphaGrid     (default [0.05 0.1 0.25 0.5 0.75 0.9 1])
%                            alpha controls lasso–ridge mixing:
%                              alpha = 1   lasso (sparse)
%                              alpha ~ 0   ridge-like (dense)
%                              0<alpha<1   elastic net
%          opts.lambdaGrid    (default logspace(-3,1,25))
%                            lambda controls overall regularization strength
%                            (larger lambda → stronger shrinkage, more zeros)
%
%        Resampling inference:
%          opts.nPerm         (default 1000) permutations for p-value
%          opts.nBoot         (default 500)  bootstrap resamples for AUC CI
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
%                             Used to compute results.AUC_global.
%
% OUTPUT (results struct)
%   Cross-validated performance (generalization estimate):
%     results.AUC               scalar   mean ROC AUC across repeats×outer folds
%     results.AUC_PR            scalar   mean precision-recall AUC across repeatsxouter folds
%     results.ACC               scalar   mean accuracy across repeats×outer folds
%     results.SENS              scalar   mean sensitivity across repeats×outer folds
%     results.SPEC              scalar   mean specificity across repeats×outer folds
%     results.ACC_balanced      scalar   mean balanced accuracy across repeatsxouter folds
%     results.allAUC            [nRepeats x outerK] fold-level ROC AUC
%     results.allAUC_PR         [nRepeats x outerK] fold-level precision-recall AUC
%     results.allACC            [nRepeats x outerK] fold-level ACC
%     results.allSENS           [nRepeats x outerK] fold-level SENS
%     results.allSPEC           [nRepeats x outerK] fold-level SPEC
%     results.allACC_balanced   [nRepeats x outerK] fold-level balanced accuracy
%
%   Model selection / weights across CV:
%     results.selectedAlpha  [nRepeats x outerK] selected alpha per outer fold
%     results.selectedLambda [nRepeats x outerK] selected lambda per outer fold
%     results.betaStore      [p x outerK x nRepeats] elastic net coefficients
%                            (feature weights only; intercept stored separately)
%     results.interceptStore [nRepeats x outerK] intercept of the selected model
%                            for each outer fold
%     results.featureWeights [p x (nRepeats*outerK)] coefficients stacked across
%                            all outer folds/repeats
%     results.meanFeatureWeight [p x 1] mean featureWeights across runs
%     results.featureStability  [p x 1] proportion of runs where |beta|>0
%
%   Global baseline (interpretation only):
%     results.AUC_global        scalar   AUC of logistic model on a global summary feature
%                                   (mean/median across features, or custom opts.globalFun)
%     results.AUC_PR_global     scalar   precision-recall AUC of logistic model on a global summary feature
%                                   (mean/median across features, or custom
%                                   opts.globalFun)
%
%   Feature stability (from CV weights):
%     results.signStability [p x 1] proportion of runs matching mean sign
%     results.selectionFrequency [p x 1] frequency of appearing in topK
%                        absolute weights across runs (topK fixed at 20)
%
%   Permutation test:
%     results.allpermAUC        [nPerm x 1] permuted ROC AUC distribution
%     results.permAUC           scalar mean permuted ROC AUC
%     results.permutation_p     scalar p = mean(permAUC >= observed AUC)
%     results.allpermAUC_PR     [nPerm x 1] permuted PR AUC distribution
%     results.permAUC_PR        scalar mean permuted PR AUC
%     results.permutation_p_PR  scalar p = mean(permAUC >= observed AUC)
%
%   Bootstrap:
%     results.allbootAUC [nBoot x 1] out-of-bag bootstrap AUC distribution
%     results.bootAUC    scalar mean out-of-bag bootstrap AUC
%     results.AUC_CI     [1 x 2] percentile CI (2.5, 97.5)
%
%   Learning curve:
%     results.learningSizes vector of sample sizes evaluated
%     results.learningAUC   vector of AUC estimates per size
%
% NOTES / INTERPRETATION (high level)
%   - Use results.AUC from nested CV as the primary generalization estimate.
%   - Even though ROC AUC is mathematically insensitive to imbalance, in
%       case of (strong) inbalance, additionally use
%       - results.AUC_PR if n(positives) is low
%       - results.balanced accuracy (mean of Sens & Spec) if n(positives) is high
%   - Elastic Net yields sparse solutions (some betas exactly zero), making:
%       featureStability (non-zero proportion) and selectionFrequency
%       particularly informative for interpretation.
%   - Bootstrap AUC is estimated using out-of-bag (OOB) testing rather than
%     evaluating on the bootstrap sample itself, which makes it more conservative
%     and typically less optimistic than naive bootstrap performance estimates.
%   - Accordingly, bootstrap AUC may be somewhat lower than nested CV AUC in
%     small samples and should be viewed primarily as a measure of sampling
%     variability rather than as the main performance estimate.
%   - Elastic net weights are sparse for alpha closer to 1 (lasso-like),
%     and dense for alpha closer to 0 (ridge-like). Stability metrics
%     (featureStability, signStability, selectionFrequency) help identify
%     robust contributors across resampling.
%
% IMPLEMENTATION NOTES
%   - Scaling is leakage-free and controlled by opts.scale.
%   - Inner CV performs a grid search over alphaGrid × lambdaGrid using AUC.
%   - Outer-fold predictions use the full elastic net linear predictor
%       score = X*beta + intercept
%     ensuring consistency with the models evaluated during inner CV.
%   - Bootstrap confidence intervals are based on out-of-bag bootstrap samples:
%     each bootstrap replicate is trained on the in-bag sample, hyperparameters
%     are tuned within that sample, and AUC is evaluated on out-of-bag subjects.
%   - If stratified cvpartition fails (e.g., extreme imbalance), the code falls back
%     to non-stratified partitions.
%   - If some folds lack a class, reduce outerK/innerK.
%
% DEPENDENCIES
%   Requires Statistics and Machine Learning Toolbox:
%     cvpartition, lassoglm, perfcurve, fitglm, grp2idx
%
% See also: lassoglm, perfcurve, cvpartition, fitglm

%% -------------------------------------------------
% 0. Defaults
%% -------------------------------------------------

if nargin < 3
    opts = struct;
end

if ~isfield(opts,'outerK'); opts.outerK = 5; end
if ~isfield(opts,'innerK'); opts.innerK = 4; end
if ~isfield(opts,'nRepeats'); opts.nRepeats = 50; end
if ~isfield(opts,'nPerm'); opts.nPerm = 1000; end
if ~isfield(opts,'nBoot'); opts.nBoot = 500; end
if ~isfield(opts,'learningSteps'); opts.learningSteps = 6; end

if ~isfield(opts,'alphaGrid')
    opts.alphaGrid = [0.05 0.1 0.25 0.5 0.75 0.9 1];
end

if ~isfield(opts,'lambdaGrid')
    opts.lambdaGrid = logspace(-3,1,25);
end

% Generic additions (kept minimal)
if ~isfield(opts,'scale'); opts.scale = 'zscore'; end       % 'zscore'|'center'|'none'
if ~isfield(opts,'globalFun'); opts.globalFun = 'mean'; end % 'mean'|'median'|function handle

rng(1,'twister')
LaBGAScore_smart_parallel_pool_setup

%% -------------------------------------------------
% 1. Outcome preparation
%% -------------------------------------------------

if iscell(Y) || isstring(Y) || iscategorical(Y)
    Y = grp2idx(Y);
end

yNum = double(Y(:)==max(Y));
[n,p] = size(X);

%% -------------------------------------------------
% 2. Repeated nested CV
%% -------------------------------------------------

AUC  = nan(opts.nRepeats,opts.outerK);
AUC_PR = nan(opts.nRepeats,opts.outerK);
ACC  = nan(opts.nRepeats,opts.outerK);
SENS = nan(opts.nRepeats,opts.outerK);
SPEC = nan(opts.nRepeats,opts.outerK);
ACC_balanced = nan(opts.nRepeats,opts.outerK);

selectedAlpha  = nan(opts.nRepeats,opts.outerK);
selectedLambda = nan(opts.nRepeats,opts.outerK);

interceptStore = nan(opts.nRepeats,opts.outerK);
betaStore = nan(p,opts.outerK,opts.nRepeats);

parfor r = 1:opts.nRepeats

    % local containers for parfor
    AUC_r  = nan(1,opts.outerK);
    AUC_PR_r  = nan(1,opts.outerK);
    ACC_r  = nan(1,opts.outerK);
    SENS_r = nan(1,opts.outerK);
    SPEC_r = nan(1,opts.outerK);
    ACC_balanced_r  = nan(1,opts.outerK);

    alpha_r  = nan(1,opts.outerK);
    lambda_r = nan(1,opts.outerK);

    intercept_r = nan(1,opts.outerK);
    beta_r = nan(p,opts.outerK);

    % robust cvpartition creation
    try
        cvOuter = cvpartition(yNum,'KFold',opts.outerK,'Stratify',true);
    catch
        cvOuter = cvpartition(n,'KFold',opts.outerK);
    end

    for k = 1:opts.outerK

        trainIdx = training(cvOuter,k);
        testIdx  = test(cvOuter,k);

        ytrain = yNum(trainIdx);
        ytest  = yNum(testIdx);

        if numel(unique(ytrain))<2 || numel(unique(ytest))<2
            continue
        end

        Xtrain = X(trainIdx,:);
        Xtest  = X(testIdx,:);

        %% leakage-free scaling (generic)
        [Xtrain, Xtest] = applyScaling(Xtrain, Xtest, opts.scale);

        %% inner CV hyperparameter search
        try
            cvInner = cvpartition(ytrain,'KFold',opts.innerK,'Stratify',true);
        catch
            cvInner = cvpartition(length(ytrain),'KFold',opts.innerK);
        end

        bestAUC = -inf;
        bestAlpha = NaN;
        bestLambda = NaN;
        bestIntercept = NaN;
        bestBeta = NaN(p,1);

        for a = 1:length(opts.alphaGrid)

            alpha = opts.alphaGrid(a);

            for f = 1:opts.innerK

                tr = training(cvInner,f);
                va = test(cvInner,f);

                ytr = ytrain(tr);
                yva = ytrain(va);

                if numel(unique(ytr))<2 || numel(unique(yva))<2
                    continue
                end

                % compute full lambda path in one call
                [B,FitInfo] = lassoglm( ...
                    Xtrain(tr,:),ytr,'binomial', ...
                    'Alpha',alpha, ...
                    'Lambda',opts.lambdaGrid, ...
                    'Standardize',false);

                for l = 1:length(opts.lambdaGrid)

                    yhat = Xtrain(va,:)*B(:,l) + FitInfo.Intercept(l);

                    [~,~,~,aucTemp] = perfcurve(yva,yhat,1);

                    if aucTemp > bestAUC

                        bestAUC = aucTemp;
                        bestAlpha = alpha;
                        bestLambda = opts.lambdaGrid(l);

                        [Bfull,FitInfoFull] = lassoglm( ...
                            Xtrain,ytrain,'binomial',...
                            'Alpha',alpha, ...
                            'Lambda',opts.lambdaGrid(l), ...
                            'Standardize',false);

                        bestBeta = Bfull;
                        bestIntercept = FitInfoFull.Intercept;

                    end

                end

            end

        end

        alpha_r(k)  = bestAlpha;
        lambda_r(k) = bestLambda;
        intercept_r(k) = bestIntercept;
        beta_r(:,k) = bestBeta;

        %% outer test prediction

        scores = Xtest*bestBeta + bestIntercept;

        [~,~,~,AUC_r(k)] = perfcurve(ytest,scores,1);

        [~, ~, ~, AUC_PR_r(k)] = perfcurve(ytest,scores,1, ...
            'xCrit','reca','yCrit','prec');

        prob = 1./(1+exp(-scores));
        pred = prob > 0.5;

        ACC_r(k) = mean(pred==ytest);

        tp = sum(pred==1 & ytest==1);
        tn = sum(pred==0 & ytest==0);
        fp = sum(pred==1 & ytest==0);
        fn = sum(pred==0 & ytest==1);

        SENS_r(k) = tp/(tp+fn);
        SPEC_r(k) = tn/(tn+fp);
        ACC_balanced_r(k) = mean([SENS_r(k),SPEC_r(k)]);

    end

    % write back to shared arrays
    AUC(r,:)  = AUC_r;
    AUC_PR(r,:)  = AUC_PR_r;
    ACC(r,:)  = ACC_r;
    SENS(r,:) = SENS_r;
    SPEC(r,:) = SPEC_r;
    ACC_balanced(r,:)  = ACC_balanced_r;

    selectedAlpha(r,:)  = alpha_r;
    selectedLambda(r,:) = lambda_r;

    interceptStore(r,:) = intercept_r;
    betaStore(:,:,r) = beta_r;

end

featureWeights = reshape(betaStore,p,opts.nRepeats*opts.outerK);

%% -------------------------------------------------
% 3. Performance summary
%% -------------------------------------------------

results.allAUC = AUC;
results.AUC = nanmean(AUC(:));

results.allAUC_PR = AUC_PR;
results.AUC_PR = nanmean(AUC_PR(:));

results.allACC = ACC;
results.ACC = nanmean(ACC(:));

results.allSENS = SENS;
results.SENS = nanmean(SENS(:));

results.allSPEC = SPEC;
results.SPEC = nanmean(SPEC(:));

results.allACC_balanced = ACC_balanced;
results.ACC_balanced = nanmean(ACC_balanced(:));

results.selectedAlpha  = selectedAlpha;
results.selectedLambda = selectedLambda;

results.interceptStore = interceptStore;
results.betaStore = betaStore;
results.featureWeights = featureWeights;

results.featureStability = mean(abs(featureWeights)>0,2);
results.meanFeatureWeight = mean(featureWeights,2);

fprintf('Nested CV AUC = %.3f\n',results.AUC)

%% -------------------------------------------------
% 4. Global signal model (generic)
%% -------------------------------------------------

if isa(opts.globalFun,'function_handle')
    globalFeature = opts.globalFun(X);
else
    switch lower(opts.globalFun)
        case 'median'
            globalFeature = median(X,2);
        otherwise
            globalFeature = mean(X,2);
    end
end

mdl = fitglm(globalFeature,yNum,'Distribution','binomial');
scores = predict(mdl,globalFeature);
[~,~,~,AUCg] = perfcurve(yNum,scores,1);
[~,~,~,AUC_PRg] = perfcurve(yNum,scores,1,...
    'xCrit','reca','yCrit','prec');
results.AUC_global = AUCg;
results.AUC_PR_global = AUC_PRg;

%% -------------------------------------------------
% 5. Sign stability
%% -------------------------------------------------

betaFlat = reshape(betaStore,p,[]);
results.signStability = mean(sign(betaFlat)==sign(results.meanFeatureWeight),2);

%% -------------------------------------------------
% 6. Top-K selection frequency
%% -------------------------------------------------

if size(X,2) < 20
    topK = size(X,2);
else
    topK = 20;
end

freq = zeros(p,1);

for i = 1:size(featureWeights,2)
    [~,idx] = sort(abs(featureWeights(:,i)),'descend');
    freq(idx(1:topK)) = freq(idx(1:topK)) + 1;
end

results.selectionFrequency = freq / size(featureWeights,2);

%% -------------------------------------------------
% 7. Permutation testing
%% -------------------------------------------------

permAUC = nan(opts.nPerm,1);

parfor i = 1:opts.nPerm
    yp = yNum(randperm(n));
    permAUC(i) = quickCV_ENet(X,yp,opts);
end

results.allpermAUC = permAUC;
results.permAUC = nanmean(permAUC);
results.permutation_p = mean(permAUC >= results.AUC,'omitnan');

figure
histogram(permAUC(~isnan(permAUC)))
hold on
xline(results.AUC)
title('Permutation AUC')

permAUC_PR = nan(opts.nPerm,1);

parfor i=1:opts.nPerm
    yp = yNum(randperm(n));
    permAUC_PR(i) = quickCV_ENet_PR(X,yp,opts);
end

results.allpermAUC_PR = permAUC_PR;
results.permAUC_PR = nanmean(permAUC_PR);
results.permutation_p_PR = mean(permAUC_PR >= results.AUC_PR,'omitnan');

figure
histogram(permAUC(~isnan(permAUC_PR)))
hold on
xline(results.AUC_PR)
title('Permutation PR AUC distribution')

%% -------------------------------------------------
% 8. Bootstrap CI (out-of-bag bootstrap)
%% -------------------------------------------------

bootAUC = nan(opts.nBoot,1);

parfor b = 1:opts.nBoot
    bootAUC(b) = bootstrapOOB_ENet(X, yNum, opts);
end

results.allbootAUC = bootAUC;
results.bootAUC = nanmean(bootAUC);
results.AUC_CI = prctile(bootAUC(~isnan(bootAUC)), [2.5 97.5]);

figure
histogram(bootAUC(~isnan(bootAUC)))
hold on
xline(results.AUC)
title('Bootstrap OOB AUC')
xlabel('AUC')
ylabel('Frequency')

%% -------------------------------------------------
% 9. Learning curve (robust stratified sampling)
%% -------------------------------------------------

sizes = round(linspace(max(10,ceil(n*0.4)), n, opts.learningSteps));
sizes = unique(sizes);

lcAUC = nan(length(sizes),1);

idxClass1 = find(yNum==1);
idxClass0 = find(yNum==0);

parfor i = 1:length(sizes)

    m = sizes(i);

    % maintain approximate class balance
    frac1 = numel(idxClass1)/n;
    n1 = max(1,round(m*frac1));
    n0 = max(1,m-n1);

    % adjust if exceeding available samples
    n1 = min(n1,length(idxClass1));
    n0 = min(n0,length(idxClass0));

    samp1 = randsample(idxClass1,n1);
    samp0 = randsample(idxClass0,n0);

    idx = [samp1; samp0];

    lcAUC(i) = quickCV_ENet(X(idx,:),yNum(idx),opts);

end

results.learningSizes = sizes;
results.learningAUC = lcAUC;

figure
plot(sizes,lcAUC,'o-')
xlabel('Sample size')
ylabel('AUC')
title('Learning curve')

end

%% -------------------------------------------------
% local helper: scaling (same as generic PLSDA)
%% -------------------------------------------------
function [XtrS, XteS] = applyScaling(Xtr, Xte, mode)

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

function AUC = bootstrapOOB_ENet(X, Y, opts)
% bootstrapOOB_ENet
% Out-of-bag bootstrap AUC for Elastic Net:
% - bootstrap sample used for training
% - OOB subjects used for testing
% - alpha/lambda selected by inner CV within the bootstrap sample
% - robustified to mirror quickCV_ENet logic

n = length(Y);

if numel(unique(Y)) < 2
    AUC = NaN;
    return
end

% -----------------------
% Bootstrap sample
% -----------------------
idxBoot = randsample(n, n, true);

% OOB subjects = never selected
inBag = false(n,1);
inBag(idxBoot) = true;
oob = ~inBag;

if sum(oob) < 2
    AUC = NaN;
    return
end

Xtrain = X(idxBoot,:);
Xtest  = X(oob,:);

ytrain = Y(idxBoot);
ytest  = Y(oob);

if numel(unique(ytrain)) < 2 || numel(unique(ytest)) < 2
    AUC = NaN;
    return
end

% -----------------------
% Leakage-free scaling
% -----------------------
[Xtrain, Xtest] = applyScaling(Xtrain, Xtest, opts.scale);

% -----------------------
% Inner CV setup
% -----------------------
innerK = min([opts.innerK, floor(numel(ytrain)/2), sum(ytrain==0), sum(ytrain==1)]);
if innerK < 2
    AUC = NaN;
    return
end

try
    cvInner = cvpartition(ytrain,'KFold',innerK,'Stratify',true);
catch
    cvInner = cvpartition(length(ytrain),'KFold',innerK);
end

% -----------------------
% Alpha grid (coarsen for speed if needed)
% -----------------------
alphaGrid = opts.alphaGrid;
if numel(alphaGrid) > 4
    alphaGrid = alphaGrid([1 round(end/2) end]);
end

bestAUC = -inf;
bestBeta = NaN(size(Xtrain,2),1);
bestIntercept = NaN;

% -----------------------
% Hyperparameter search
% -----------------------
for a = 1:numel(alphaGrid)

    alpha = alphaGrid(a);

    % Fit lambda path once for this alpha on full bootstrap training set
    try
        [BfullPath,FitInfoFullPath] = lassoglm( ...
            Xtrain, ytrain, 'binomial', ...
            'Alpha', alpha, ...
            'Standardize', false, ...
            'NumLambda', 25, ...
            'LambdaRatio', 1e-3, ...
            'MaxIter', 1e4, ...
            'RelTol', 1e-3);
    catch
        continue
    end

    lambdaPath = FitInfoFullPath.Lambda;
    nLambda = numel(lambdaPath);

    foldAUC = nan(innerK, nLambda);

    for f = 1:innerK

        tr = training(cvInner,f);
        va = test(cvInner,f);

        ytr = ytrain(tr);
        yva = ytrain(va);

        if numel(unique(ytr))<2 || numel(unique(yva))<2
            continue
        end

        % Recompute path on inner-training fold using same lambda path
        try
            [B,FitInfo] = lassoglm( ...
                Xtrain(tr,:), ytr, 'binomial', ...
                'Alpha', alpha, ...
                'Lambda', lambdaPath, ...
                'Standardize', false, ...
                'MaxIter', 1e4, ...
                'RelTol', 1e-3);
        catch
            continue
        end

        for l = 1:nLambda
            score = Xtrain(va,:)*B(:,l) + FitInfo.Intercept(l);

            try
                [~,~,~,foldAUC(f,l)] = perfcurve(yva, score, 1);
                if ~isfinite(foldAUC(f,l))
                    foldAUC(f,l) = NaN;
                end
            catch
                foldAUC(f,l) = NaN;
            end
        end
    end

    meanAUC = nanmean(foldAUC,1);

    if all(isnan(meanAUC))
        continue
    end

    % -----------------------
    % Robust lambda selection
    % -----------------------
    % First choose min-deviance over inner folds
    [~,idxm] = max(meanAUC);

    % Approximate 1SE-style choice:
    % pick the most regularized lambda whose AUC is within 1 SE of the best
    seAUC = nanstd(foldAUC,[],1) ./ sqrt(sum(~isnan(foldAUC),1));
    bestMean = meanAUC(idxm);
    thresh = bestMean - seAUC(idxm);

    idxCandidates = find(meanAUC >= thresh);

    if ~isempty(idxCandidates)
        idx1 = idxCandidates(1);   % most regularized among acceptable
    else
        idx1 = idxm;
    end

    idx = idx1;

    % If 1SE choice is all-zero, try min-deviance
    if all(BfullPath(:,idx)==0) && ~isempty(idxm) && idxm>=1
        idx = idxm;
    end

    % If still all-zero, choose best non-zero solution on path
    if all(BfullPath(:,idx)==0)
        nonzeroCols = find(any(BfullPath~=0,1));
        validCols = nonzeroCols(isfinite(meanAUC(nonzeroCols)));
        if ~isempty(validCols)
            [~,ii] = max(meanAUC(validCols));
            idx = validCols(ii);
        end
    end

    if ~isfinite(meanAUC(idx))
        continue
    end

    if meanAUC(idx) > bestAUC
        bestAUC = meanAUC(idx);
        bestBeta = BfullPath(:,idx);
        bestIntercept = FitInfoFullPath.Intercept(idx);
    end

end

% -----------------------
% Final OOB evaluation
% -----------------------
if all(isnan(bestBeta)) || isnan(bestIntercept)
    AUC = NaN;
    return
end

score = Xtest*bestBeta + bestIntercept;

try
    [~,~,~,AUC] = perfcurve(ytest, score, 1);
    if ~isfinite(AUC)
        AUC = NaN;
    end
catch
    AUC = NaN;
end

end