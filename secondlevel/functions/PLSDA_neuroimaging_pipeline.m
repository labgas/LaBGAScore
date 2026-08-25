function results = PLSDA_neuroimaging_pipeline(X,Y,opts)

% PLSDA_neuroimaging_pipeline  Robust PLS-DA pipeline for neuroimaging feature matrices.
%
% This function implements Partial Least Squares Discriminant Analysis (PLS-DA)
% for binary classification with repeated nested k-fold cross-validation.
% It is designed for neuroimaging feature matrices (subjects × features) such
% as PET ROI binding, fMRI ROI betas, morphometry, connectivity-derived
% measures (edges), graph metrics, or multimodal ROI feature sets. The architecture emphasizes:
%   - leakage-free preprocessing (train-only scaling in every fold)
%   - inner CV selection of number of latent variables (LVs)
%   - outer CV estimation of generalization performance
%   - resampling-based inference (permutation p-value, bootstrap CI)
%   - stability/interpretability metrics (VIP, stabilityZ, sign stability,
%     top-K selection frequency)
%   - learning curve (AUC vs sample size)
%
% USAGE
%   results = PLSDA_neuroimaging_pipeline(X, Y)
%   results = PLSDA_neuroimaging_pipeline(X, Y, opts)
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
%                             Used for results.AUC_global and the *_global_cv fields.
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
%
%                             NOTE for classification: Y stays binary, so only X
%                             is residualized. If a covariate is itself
%                             associated with the class, removing it from X also
%                             removes part of the class signal. That is
%                             conservative by design, not a defect.
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
%     results.selectedLV [nRepeats x outerK] selected LV per outer fold
%     results.betaStore  [p+1 x outerK x nRepeats] PLS regression betas
%                        (row 1 is intercept; rows 2..p+1 correspond to features)
%     results.featureWeights [p x (nRepeats*outerK)] betas (no intercept),
%                        stacked across all outer folds/repeats
%     results.meanFeatureWeight [p x 1] mean featureWeights across runs
%
%   Global baseline:
%     results.AUC_global        scalar   IN-SAMPLE ROC AUC of a logistic model on
%                                   a global summary feature, fit and predicted on
%                                   all subjects. Retained unchanged for
%                                   continuity, but NOT comparable to the
%                                   cross-validated results.AUC.
%     results.AUC_global_cv     scalar   cross-validated counterpart, through the
%                                   same repeated outer CV as the model. THIS is
%                                   what to compare against results.AUC.
%     results.AUC_PR_global_cv  scalar   cross-validated PR-AUC counterpart.
%     results.AUC_PR_global     scalar   precision-recall AUC of logistic model on a global summary feature
%                                   (mean/median across features, or custom
%                                   opts.globalFun)
%
%   Final model on all data (interpretation only; NOT for performance):
%     results.finalLV        scalar median selected LV across all runs (then capped valid)
%     results.betaFinal      [p+1 x 1] betas from final model
%     results.varExplainedX  [1 x finalLV] PCTVAR(1,:) from plsregress
%     results.varExplainedY  [1 x finalLV] PCTVAR(2,:) from plsregress
%     results.finalXLoadings [p x finalLV] XL from plsregress
%     results.finalYLoadings [1 x finalLV] YL from plsregress
%
%   Feature importance / stability (from final model and CV betas):
%     results.VIP         [p x 1] VIP scores (higher = more important)
%     results.meanBeta    [p x 1] mean beta across outer CV runs
%     results.sdBeta      [p x 1] std beta across outer CV runs
%     results.stabilityZ  [p x 1] meanBeta ./ sdBeta (stability statistic)
%     results.signStability [p x 1] proportion of runs matching mean sign
%
%   Permutation test (label permutation):
%     results.allpermAUC        [nPerm x 1] permuted ROC AUC distribution
%     results.permAUC           scalar mean permuted ROC AUC
%     results.quickCV_observed  scalar observed AUC from quickCV_PLSDA, i.e. the SAME
%                            estimator that generates the null
%     results.permutation_p     scalar (sum(permAUC >= quickCV_observed) + 1) / (nValid + 1)
%     results.allpermAUC_PR     [nPerm x 1] permuted PR AUC distribution
%     results.permAUC_PR        scalar mean permuted PR AUC
%     results.quickCV_observed_PR  scalar matched observed PR AUC
%     results.permutation_p_PR  scalar as above, for PR AUC
%
%     Labels are permuted directly. When covariates are supplied, X is
%     residualized on them inside every training fold, so the features the model
%     sees are already orthogonal to the covariates and free label permutation
%     is a valid test of the partial X-Y association (the Kennedy scheme).
%
%     Both metrics come from the SAME permutation draw, so the two p-values are
%     mutually coherent and the stage costs one pass rather than two.
%
%     The p-value uses the matched observed statistic rather than results.AUC
%     because the null comes from quickCV_PLSDA. Mixing a repeated-nested-CV AUC with a
%     quickCV null gave a measured false positive rate near 40 percent at
%     alpha = 0.05. results.AUC remains the performance estimate to report.
%
%   Bootstrap:
%     results.allbootAUC [nBoot x 1] out-of-bag bootstrap ROC AUC distribution
%     results.bootAUC    scalar mean out-of-bag bootstrap ROC AUC
%     results.AUC_CI     [1 x 2] percentile CI (2.5, 97.5)
%
%   Learning curve:
%     results.learningSizes vector of sample sizes evaluated
%     results.learningAUC   vector of AUC estimates per size
%
%   Provenance:
%     results.covariateInfo struct: .used .names .nCov .rank .residualizeY
%                        .order ('residualize-then-scale') .permScheme
%     results.seed          the RNG seed actually used
%
% NOTES / INTERPRETATION (high level)
%   - Use results.AUC from nested CV as the primary generalization estimate.
%   - Even though ROC AUC is mathematically insensitive to imbalance, in
%       case of (strong) inbalance, additionally use
%       - results.AUC_PR if n(positives) is low
%       - results.balanced accuracy (mean of Sens & Spec) if n(positives) is high 
%   - Bootstrap AUC is estimated using out-of-bag (OOB) testing rather than
%     evaluating on the bootstrap sample itself, which makes it more conservative
%     and typically less optimistic than naive bootstrap performance estimates.
%   - Accordingly, bootstrap AUC may be somewhat lower than nested CV AUC in
%     small samples and should be viewed primarily as a measure of sampling
%     variability rather than as the main performance estimate.
%   - VIP and stabilityZ provide complementary interpretability:
%       VIP: importance in the final fitted model
%       stabilityZ: robustness of the effect across CV runs
%
% IMPLEMENTATION NOTES
%   - Scaling is leakage-free and controlled by opts.scale.
%   - Bootstrap confidence intervals are based on out-of-bag bootstrap samples:
%     each bootstrap replicate is trained on the in-bag sample, latent variables
%     are tuned within that sample, and AUC is evaluated on out-of-bag subjects.
%   - If stratified cvpartition fails (e.g., extreme imbalance), the code falls back
%     to non-stratified partitions.
%   - If some folds lack a class, reduce outerK/innerK.
%
% DEPENDENCIES
%   Requires Statistics and Machine Learning Toolbox:
%     cvpartition, plsregress, perfcurve, fitglm, grp2idx
%
% See also: plsregress, perfcurve, cvpartition, fitglm

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
if ~isfield(opts,'seed'); opts.seed = 1; end                % see setParforStream

warnUnknownOptions(opts, { ...
    'outerK','innerK','nRepeats','maxLV','nPerm','nBoot','learningSteps', ...
    'scale','globalFun','seed','covariates','covariateNames'}, ...
    'PLSDA_neuroimaging_pipeline');

rng(opts.seed,'twister')

%% -------------------------------------------------
% 1. Outcome preparation
%% -------------------------------------------------

if iscell(Y) || isstring(Y) || iscategorical(Y)
    Y = grp2idx(Y);
end

yNum = double(Y(:)==max(Y));
[n,p] = size(X);

if numel(yNum) ~= n
    error('PLSDA_neuroimaging_pipeline:sizeMismatch', ...
        'X has %d rows but Y has %d elements.', n, numel(yNum));
end

if ~all(isfinite(X(:)))
    error('PLSDA_neuroimaging_pipeline:nonFiniteX', ...
        'X contains NaN or Inf. Handle missing values upstream (impute or drop).');
end

% Covariates are resolved ONCE here and then passed explicitly to every helper
% alongside X; see foldPreprocess for why they never travel inside opts.
[Cov, covInfo] = validateCovariates(opts, n, 'PLSDA_neuroimaging_pipeline');

%% -------------------------------------------------
% 2. Repeated Nested Cross-Validation
%% -------------------------------------------------

AUC  = nan(opts.nRepeats,opts.outerK);
AUC_PR = nan(opts.nRepeats,opts.outerK);
ACC  = nan(opts.nRepeats,opts.outerK);
SENS = nan(opts.nRepeats,opts.outerK);
SPEC = nan(opts.nRepeats,opts.outerK);
ACC_balanced = nan(opts.nRepeats,opts.outerK);

selectedLV = nan(opts.nRepeats,opts.outerK);
betaStore = nan(p+1,opts.outerK,opts.nRepeats);
featureWeights = nan(p, opts.nRepeats*opts.outerK);

for r = 1:opts.nRepeats

    % robust cvpartition creation
    try
        cvOuter = cvpartition(yNum,'KFold',opts.outerK,'Stratify',true);
    catch
        cvOuter = cvpartition(length(yNum),'KFold',opts.outerK);
    end

    for k = 1:opts.outerK

        trainIdx = training(cvOuter,k);
        testIdx  = test(cvOuter,k);

        ytrain = yNum(trainIdx);
        ytest  = yNum(testIdx);

        % skip if class missing
        if numel(unique(ytrain))<2 || numel(unique(ytest))<2
            continue
        end

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

        %% inner CV LV tuning

        maxLV = capLV(opts.maxLV, Xtrain);

        if maxLV < 1
            continue
        end

        innerAUC = nan(maxLV,1);

        try
            cvInner = cvpartition(ytrain,'KFold',opts.innerK,'Stratify',true);
        catch
            cvInner = cvpartition(length(ytrain),'KFold',opts.innerK);
        end

        parfor lv = 1:maxLV

            foldAUC = nan(opts.innerK,1);

            for f = 1:opts.innerK

                tr = training(cvInner,f);
                va = test(cvInner,f);

                ytr = ytrain(tr);
                yva = ytrain(va);

                if numel(unique(ytr))<2 || numel(unique(yva))<2
                    continue
                end

                % Each inner fold preprocesses itself from the RAW training
                % rows; reusing the outer fold's constants would leak
                % outer-training statistics (and the outer nuisance fit) into
                % every inner validation fold.
                C2 = []; Cv = [];
                if ~isempty(Ctrain), C2 = Ctrain(tr,:); Cv = Ctrain(va,:); end

                [X2, Xv] = foldPreprocess(Xtrain_raw(tr,:), Xtrain_raw(va,:), C2, Cv, opts.scale);

                if capLV(lv, X2) < lv
                    continue
                end

                [~,~,~,~,beta] = plsregress(X2,ytr,lv);

                yhat = [ones(sum(va),1) Xv]*beta;

                [~,~,~,foldAUC(f)] = perfcurve(yva,yhat,1);

            end

            innerAUC(lv) = mean(foldAUC,'omitnan');

        end

        if all(isnan(innerAUC))
            continue
        end
        [~,bestLV] = max(innerAUC);
        selectedLV(r,k) = bestLV;

        %% fit model

        [~,~,~,~,beta] = plsregress(Xtrain,ytrain,bestLV);
        betaStore(:,k,r) = beta;
        featureWeights(:,(r-1)*opts.outerK+k) = beta(2:end);

        scores = [ones(sum(testIdx),1) Xtest]*beta;

        [~,~,~,AUC(r,k)] = perfcurve(ytest,scores,1);

        [~, ~, ~, AUC_PR(r,k)] = perfcurve(ytest,scores,1, ...
            'xCrit','reca','yCrit','prec');


        pred = scores>0.5;

        ACC(r,k) = mean(pred==ytest);

        tp = sum(pred==1 & ytest==1);
        tn = sum(pred==0 & ytest==0);
        fp = sum(pred==1 & ytest==0);
        fn = sum(pred==0 & ytest==1);

        SENS(r,k) = tp/(tp+fn);
        SPEC(r,k) = tn/(tn+fp);
        ACC_balanced(r,k) = mean([SENS(r,k),SPEC(r,k)]);

    end
end

results.allAUC  = AUC;
results.AUC     = mean(AUC(:),'omitnan');
results.allAUC_PR  = AUC_PR;
results.AUC_PR     = mean(AUC_PR(:),'omitnan');
results.allACC  = ACC;
results.ACC     = mean(ACC(:),'omitnan');
results.allSENS = SENS;
results.SENS    = mean(SENS(:),'omitnan');
results.allSPEC = SPEC;
results.SPEC    = mean(SPEC(:),'omitnan');
results.allACC_balanced = ACC_balanced;
results.ACC_balanced = mean(ACC_balanced(:),'omitnan');

results.selectedLV    = selectedLV;
results.betaStore     = betaStore;
results.featureWeights = featureWeights;

% Skipped folds (missing class) leave all-NaN columns; omit them rather than
% letting a single skipped fold turn the whole weight vector into NaN.
results.meanFeatureWeight  = mean(featureWeights,2,'omitnan');

fprintf('Nested CV AUC = %.3f\n',results.AUC)

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
            globalFeature = mean(X,2); % default (matches your TSPO script)
    end
end

mdl = fitglm(globalFeature,yNum,'Distribution','binomial');
scores = predict(mdl,globalFeature);
[~,~,~,AUCg] = perfcurve(yNum,scores,1);
[~,~,~,AUC_PRg] = perfcurve(yNum,scores,1,...
    'xCrit','reca','yCrit','prec');
results.AUC_global = AUCg;
results.AUC_PR_global = AUC_PRg;

% Cross-validated companions. AUC_global / AUC_PR_global above are IN-SAMPLE
% (fitglm predicts the same subjects it was fit on) and so are not comparable
% to the cross-validated results.AUC; the *_cv fields are.
gcv = globalBaselineCV(X, yNum, Cov, opts, 'classify');
results.AUC_global_cv    = gcv.AUC;
results.AUC_PR_global_cv = gcv.AUC_PR;

%% -------------------------------------------------
% 4. Final model (interpretation only)
%% -------------------------------------------------

% apply same scaling on all data (generic)
% Interpretation-only fit on all subjects; no held-out set, so a full-sample
% nuisance fit is correct by construction.
[Xz,~] = foldPreprocess(X, [], Cov, [], opts.scale);

finalLV = round(median(selectedLV(:),'omitnan'));
finalLV = max(1, min(finalLV, capLV(opts.maxLV, Xz)));

[XL,YL,XS,~,beta,PCTVAR,~,stats] = plsregress(Xz,yNum,finalLV);

results.finalLV = finalLV;
results.betaFinal = beta;
results.varExplainedX = PCTVAR(1,:);
results.varExplainedY = PCTVAR(2,:);
results.finalXLoadings = XL;
results.finalYLoadings = YL;

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

% The null comes from quickCV_PLSDA, so the observed statistic it is compared
% against must come from quickCV_PLSDA too. Comparing the repeated-nested-CV
% AUC against a quickCV null mixes two different estimators and, measured on
% true-null data, produced a false positive rate near 40% at alpha = 0.05.
%
% Labels are permuted directly. With covariates supplied, X is residualized on
% them inside every training fold, so the features the model sees are already
% orthogonal to the covariates and free label permutation is a valid test of
% the partial X-Y association (the Kennedy scheme). Without covariates this is
% the usual unrestricted permutation.
obsMatched    = quickCV_PLSDA(X, yNum, Cov, opts);
obsMatched_PR = quickCV_PLSDA_PR(X, yNum, Cov, opts);

results.quickCV_observed    = obsMatched;
results.quickCV_observed_PR = obsMatched_PR;

permAUC    = nan(opts.nPerm,1);
permAUC_PR = nan(opts.nPerm,1);

% One permutation draw feeds both metrics: the ROC and PR nulls are then
% mutually coherent, and this halves the permutation runtime.
parfor i=1:opts.nPerm
    setParforStream(opts.seed, 2e6, i);
    yp = yNum(randperm(n));
    permAUC(i)    = quickCV_PLSDA(X, yp, Cov, opts);
    permAUC_PR(i) = quickCV_PLSDA_PR(X, yp, Cov, opts);
end

results.allpermAUC = permAUC;
results.permAUC = mean(permAUC,'omitnan');
results.permutation_p = (sum(permAUC >= obsMatched) + 1) / (sum(~isnan(permAUC)) + 1);

figure
histogram(permAUC(~isnan(permAUC)))
hold on
xline(obsMatched)
title('Permutation ROC AUC distribution')

results.allpermAUC_PR = permAUC_PR;
results.permAUC_PR = mean(permAUC_PR,'omitnan');
results.permutation_p_PR = (sum(permAUC_PR >= obsMatched_PR) + 1) / (sum(~isnan(permAUC_PR)) + 1);

figure
histogram(permAUC_PR(~isnan(permAUC_PR)))
hold on
xline(obsMatched_PR)
title('Permutation PR AUC distribution')

%% -------------------------------------------------
% 8. Bootstrap AUC CI (out-of-bag bootstrap)
%% -------------------------------------------------

bootAUC = nan(opts.nBoot,1);

parfor b = 1:opts.nBoot
    setParforStream(opts.seed, 3e6, b);
    bootAUC(b) = bootstrapOOB_PLSDA(X, yNum, Cov, opts);
end

results.allbootAUC = bootAUC;
results.bootAUC = mean(bootAUC,'omitnan');
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

    setParforStream(opts.seed, 4e6, i);

    m = sizes(i);

    frac1 = numel(idxClass1)/n;
    n1 = max(1,round(m*frac1));
    n0 = max(1,m-n1);

    n1 = min(n1,length(idxClass1));
    n0 = min(n0,length(idxClass0));

    samp1 = randsample(idxClass1,n1);
    samp0 = randsample(idxClass0,n0);

    idx = [samp1; samp0];

    Cidx = []; if ~isempty(Cov), Cidx = Cov(idx,:); end
    lcAUC(i) = quickCV_PLSDA(X(idx,:),yNum(idx),Cidx,opts);

end

results.learningSizes = sizes;
results.learningAUC = lcAUC;

figure
plot(sizes,lcAUC,'o-')
xlabel('Sample size')
ylabel('AUC')
title('Learning curve')

%% ---------------------------
% 10. Sign Stability
% ----------------------------

% betaStore is (p+1) x outerK x nRepeats, so compare sign across folds/repeats
betaNoIntercept = betaStore(2:end,:,:);
betaFlat = reshape(betaNoIntercept,p,[]);
% Compare each fold/repeat sign to the sign of the mean weight
signStability = mean(sign(results.meanFeatureWeight) == sign(betaFlat),2,'omitnan');
results.signStability = signStability;

%% -------------------------------------------------
% Provenance
%% -------------------------------------------------

covInfo.residualizeY = false;   % never applicable: Y is binary here
covInfo.order        = 'residualize-then-scale';
covInfo.permScheme   = 'label-permutation';
results.covariateInfo = covInfo;
results.seed          = opts.seed;

end
