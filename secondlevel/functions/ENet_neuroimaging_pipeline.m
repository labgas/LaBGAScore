function results = ENet_neuroimaging_pipeline(X,Y,opts)
% ENet_neuroimaging_pipeline  Robust Elastic Net pipeline for neuroimaging feature matrices.
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
%        Elastic Net specifics:
%          opts.tuneRule      (default '1se') inner-CV selection rule:
%                             '1se' picks the most regularized (alpha, lambda)
%                                   whose fold-averaged AUC is within one
%                                   standard error of the best -- sparser and
%                                   more stable weight maps;
%                             'max' picks the plain fold-averaged argmax.
%                             Measured on data with real signal, '1se' cost
%                             about 0.01 AUC against 'max' while selecting far
%                             stronger regularization.
%          opts.selectionTopK (default min(20, max(3, ceil(0.25*p))), capped at
%                             p-1) size of the top-K set for selectionFrequency.
%          opts.doFinalModel  (default false) additionally fit an
%                             interpretation-only Elastic Net on all subjects at
%                             the median selected (alpha, lambda), adding
%                             results.betaFinal / .interceptFinal / .finalAlpha
%                             / .finalLambda. Off by default because for a
%                             sparse model the mean of the CV betas is a
%                             defensible bagged weight map, whereas a single
%                             full-data fit commits to one arbitrary support set.
%          opts.verbose       (default false) print per-fold diagnostics from
%                             quickCV_ENet. Leave false: with parfor these
%                             stream from every worker on every permutation.
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
%   Feature stability (from CV weights):
%     results.signStability [p x 1] proportion of runs matching mean sign
%     results.selectionFrequency [p x 1] frequency of appearing in the topK
%                        absolute weights across valid runs. topK defaults to
%                        min(20, max(3, ceil(0.25*p))) capped at p-1, NOT a
%                        fixed 20; override with opts.selectionTopK. Reported
%                        back in results.selectionTopK.
%     results.nValidFolds scalar number of outer folds that produced a model
%                        (skipped folds are excluded from the stability metrics
%                        rather than counted as "feature not selected")
%
%   Permutation test (label permutation):
%     results.allpermAUC        [nPerm x 1] permuted ROC AUC distribution
%     results.permAUC           scalar mean permuted ROC AUC
%     results.quickCV_observed  scalar observed AUC from quickCV_ENet, i.e. the SAME
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
%     because the null comes from quickCV_ENet. Mixing a repeated-nested-CV AUC with a
%     quickCV null gave a measured false positive rate near 40 percent at
%     alpha = 0.05. results.AUC remains the performance estimate to report.
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

if ~isfield(opts,'selectionTopK') || isempty(opts.selectionTopK)
    p = size(X,2);
    topK = min(20, max(3, ceil(0.25 * p)));
    topK = min(topK, p-1);
else
    topK = min(opts.selectionTopK, size(X,2)-1);
    topK = max(1, topK);
end

% Generic additions (kept minimal)
if ~isfield(opts,'scale'); opts.scale = 'zscore'; end       % 'zscore'|'center'|'none'
if ~isfield(opts,'globalFun'); opts.globalFun = 'mean'; end % 'mean'|'median'|function handle
if ~isfield(opts,'tuneRule'); opts.tuneRule = '1se'; end    % '1se'|'max', see selectENetHyperparams
if ~isfield(opts,'doFinalModel'); opts.doFinalModel = false; end % optional full-data interpretive fit
if ~isfield(opts,'verbose'); opts.verbose = false; end      % per-fold diagnostics from quickCV_ENet
if ~isfield(opts,'seed'); opts.seed = 1; end                % RNG seed, see setParforStream

warnUnknownOptions(opts, { ...
    'outerK','innerK','nRepeats','nPerm','nBoot','learningSteps', ...
    'alphaGrid','lambdaGrid','selectionTopK','scale','globalFun', ...
    'tuneRule','verbose','doFinalModel','seed','covariates','covariateNames'}, ...
    'ENet_neuroimaging_pipeline');

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
    error('ENet_neuroimaging_pipeline:sizeMismatch', ...
        'X has %d rows but Y has %d elements.', n, numel(yNum));
end

if ~all(isfinite(X(:)))
    error('ENet_neuroimaging_pipeline:nonFiniteX', ...
        'X contains NaN or Inf. Handle missing values upstream (impute or drop).');
end

% Covariates are resolved ONCE here and then passed explicitly to every helper
% alongside X. They are deliberately never re-read from opts inside a helper:
% the bootstrap and learning-curve stages subset subjects, and a full-length
% covariate matrix reaching a subsetted X would misalign silently.
[Cov, covInfo] = validateCovariates(opts, n, 'ENet_neuroimaging_pipeline');

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

    setParforStream(opts.seed, 1e6, r);

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

        %% inner CV hyperparameter search
        try
            cvInner = cvpartition(ytrain,'KFold',opts.innerK,'Stratify',true);
        catch
            cvInner = cvpartition(length(ytrain),'KFold',opts.innerK);
        end

        % Hyperparameters are selected on the FOLD-AVERAGED inner AUC, and the
        % model is refit exactly once at the winning (alpha, lambda). The former
        % inline search compared a running maximum against a single inner fold's
        % AUC over a flattened alpha x fold x lambda grid, so the winner was the
        % luckiest fold rather than the best-performing pair.
        % Raw training features go in: each inner fold preprocesses itself.
        [bestBeta, bestIntercept, bestAlpha, bestLambda] = selectENetHyperparams( ...
            Xtrain_raw, ytrain, Ctrain, cvInner, opts.alphaGrid, ...
            enetLambdaGrid(opts), opts.tuneRule, opts.scale);

        % No inner fold produced a usable model: leave this outer fold as NaN
        % rather than scoring with a NaN beta (which errors inside perfcurve).
        if ~isfinite(bestIntercept) || all(isnan(bestBeta))
            continue
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
results.AUC = mean(AUC(:),'omitnan');

results.allAUC_PR = AUC_PR;
results.AUC_PR = mean(AUC_PR(:),'omitnan');

results.allACC = ACC;
results.ACC = mean(ACC(:),'omitnan');

results.allSENS = SENS;
results.SENS = mean(SENS(:),'omitnan');

results.allSPEC = SPEC;
results.SPEC = mean(SPEC(:),'omitnan');

results.allACC_balanced = ACC_balanced;
results.ACC_balanced = mean(ACC_balanced(:),'omitnan');

results.selectedAlpha  = selectedAlpha;
results.selectedLambda = selectedLambda;

results.interceptStore = interceptStore;
results.betaStore = betaStore;
results.featureWeights = featureWeights;

% Folds that were skipped (missing class, degenerate fit) leave an all-NaN
% column in featureWeights. Mask those out rather than letting them count as
% "feature not selected" (NaN>0 is false) or poison the mean.
validFolds = ~all(isnan(featureWeights),1);

if any(validFolds)
    results.featureStability  = mean(abs(featureWeights(:,validFolds))>0,2,'omitnan');
    results.meanFeatureWeight = mean(featureWeights(:,validFolds),2,'omitnan');
else
    results.featureStability  = nan(p,1);
    results.meanFeatureWeight = nan(p,1);
end

results.nValidFolds = sum(validFolds);

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

% Cross-validated companions; see PLSDA_neuroimaging_pipeline for the rationale.
gcv = globalBaselineCV(X, yNum, Cov, opts, 'classify');
results.AUC_global_cv    = gcv.AUC;
results.AUC_PR_global_cv = gcv.AUC_PR;

%% -------------------------------------------------
% 5. Sign stability
%% -------------------------------------------------

betaFlat = reshape(betaStore,p,[]);
betaFlat = betaFlat(:,~all(isnan(betaFlat),1));
results.signStability = mean(sign(betaFlat)==sign(results.meanFeatureWeight),2,'omitnan');

%% -------------------------------------------------
% 6. Top-K selection frequency
%% -------------------------------------------------

freq = zeros(size(results.featureWeights,1),1);

% Only count folds that actually produced weights; sorting an all-NaN column
% would otherwise increment an arbitrary set of indices.
freqCols = find(validFolds);

for i = freqCols(:)'
    [~,idx] = sort(abs(results.featureWeights(:,i)),'descend','MissingPlacement','last');
    freq(idx(1:topK)) = freq(idx(1:topK)) + 1;
end

freq = freq / max(numel(freqCols),1);
results.selectionFrequency = freq;
results.selectionTopK = topK;

%% -------------------------------------------------
% 6b. Final full-data model (interpretation only, optional)
%% -------------------------------------------------

% The PLS pipelines all fit an interpretation-only model on the full sample;
% ENet historically did not, so plot_ENet_diagnostics_neuroimaging works from
% the CV-averaged weights instead. That is arguably the better choice for a
% sparse model -- the mean of CV betas is a legitimate bagged weight map, while
% a single full-data fit commits to one arbitrary support set. The option is
% therefore OFF by default, which keeps the results struct unchanged.
if opts.doFinalModel

    [Xz,~] = foldPreprocess(X, [], Cov, [], opts.scale);

    finalAlpha  = median(selectedAlpha(:),  'omitnan');
    finalLambda = median(selectedLambda(:), 'omitnan');

    if isfinite(finalAlpha) && isfinite(finalLambda)
        % snap alpha back onto the supplied grid; the median of a discrete grid
        % need not itself be a grid value
        [~,ia] = min(abs(opts.alphaGrid - finalAlpha));
        finalAlpha = opts.alphaGrid(ia);

        try
            [Bf, FIf] = lassoglm(Xz, yNum, 'binomial', ...
                'Alpha', finalAlpha, 'Lambda', finalLambda, 'Standardize', false);
            results.betaFinal      = Bf(:,1);
            results.interceptFinal = FIf.Intercept(1);
            results.finalAlpha     = finalAlpha;
            results.finalLambda    = finalLambda;
        catch ME
            warning('ENet_neuroimaging_pipeline:finalModelFailed', ...
                'Final full-data fit failed (%s); results.betaFinal not set.', ME.message);
        end
    end

end

%% -------------------------------------------------
% 7. Permutation testing
%% -------------------------------------------------

% The p-value compares like with like. The null is generated by quickCV_ENet,
% so the observed statistic it is compared against must come from quickCV_ENet
% too. Comparing the repeated-nested-CV AUC against a quickCV null mixes two
% different estimators and, measured on true-null data, produced a false
% positive rate near 40% at alpha = 0.05.
%
% Labels are permuted directly. With covariates supplied, X is residualized on
% them inside every training fold, so the features the model sees are already
% orthogonal to the covariates and free label permutation is a valid test of
% the partial X-Y association (the Kennedy scheme). Without covariates this is
% simply the usual unrestricted permutation.
obsMatched    = quickCV_ENet(X, yNum, Cov, opts);
obsMatched_PR = quickCV_ENet_PR(X, yNum, Cov, opts);

results.quickCV_observed    = obsMatched;
results.quickCV_observed_PR = obsMatched_PR;

permAUC    = nan(opts.nPerm,1);
permAUC_PR = nan(opts.nPerm,1);

% One permutation draw feeds both metrics: the ROC and PR nulls are then
% mutually coherent, and this halves the permutation runtime.
parfor i = 1:opts.nPerm
    setParforStream(opts.seed, 2e6, i);
    yp = yNum(randperm(n));
    permAUC(i)    = quickCV_ENet(X, yp, Cov, opts);
    permAUC_PR(i) = quickCV_ENet_PR(X, yp, Cov, opts);
end

results.allpermAUC = permAUC;
results.permAUC = mean(permAUC,'omitnan');
results.permutation_p = (sum(permAUC >= obsMatched) + 1) / (sum(~isnan(permAUC)) + 1);

figure
histogram(permAUC(~isnan(permAUC)))
hold on
xline(obsMatched)
title('Permutation AUC')

results.allpermAUC_PR = permAUC_PR;
results.permAUC_PR = mean(permAUC_PR,'omitnan');
results.permutation_p_PR = (sum(permAUC_PR >= obsMatched_PR) + 1) / (sum(~isnan(permAUC_PR)) + 1);

figure
histogram(permAUC_PR(~isnan(permAUC_PR)))
hold on
xline(obsMatched_PR)
title('Permutation PR AUC distribution')

%% -------------------------------------------------
% 8. Bootstrap CI (out-of-bag bootstrap)
%% -------------------------------------------------

bootAUC = nan(opts.nBoot,1);

parfor b = 1:opts.nBoot
    setParforStream(opts.seed, 3e6, b);
    bootAUC(b) = bootstrapOOB_ENet(X, yNum, Cov, opts);
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

    Cidx = []; if ~isempty(Cov), Cidx = Cov(idx,:); end
    lcAUC(i) = quickCV_ENet(X(idx,:),yNum(idx),Cidx,opts);

end

results.learningSizes = sizes;
results.learningAUC = lcAUC;

figure
plot(sizes,lcAUC,'o-')
xlabel('Sample size')
ylabel('AUC')
title('Learning curve')

%% -------------------------------------------------
% 10. Provenance
%% -------------------------------------------------

covInfo.residualizeY = false;   % never applicable: Y is binary here
covInfo.order        = 'residualize-then-scale';
covInfo.permScheme   = 'label-permutation';
results.covariateInfo = covInfo;
results.tuneRule      = opts.tuneRule;
results.seed          = opts.seed;

end
