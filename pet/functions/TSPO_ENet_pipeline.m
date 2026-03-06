function results = TSPO_ENet_pipeline(X,Y,opts)

% Robust Elastic Net pipeline for neuroimaging ROI features.
%
% This function implements Elastic Net regularized logistic regression for
% binary classification with repeated nested k-fold cross-validation.
% It is designed for neuroimaging feature matrices (subjects × features) such
% as PET ROI binding, fMRI ROI betas, morphometry, or connectivity-derived
% measures. The architecture emphasizes:
%   - leakage-free preprocessing (train-only scaling in every fold)
%   - inner CV selection of Elastic Net hyperparameters (alpha, lambda)
%   - outer CV estimation of generalization performance
%   - resampling-based inference (permutation p-value, bootstrap CI)
%   - stability/interpretability metrics (sparsity/stability, sign stability,
%     top-K selection frequency)
%   - learning curve (AUC vs sample size)
%
% USAGE
%   results = TSPO_ENet_pipeline(X, Y)
%   results = TSPO_ENet_pipeline(X, Y, opts)
%
% INPUTS
%   X    [n x p] numeric
%        Feature matrix (n subjects, p features). Each column is a feature
%        (e.g., ROI binding, ROI beta, edge weight, graph metric). Missing
%        values should be handled upstream (impute/remove) prior to calling
%        this function.
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
%                            (larger lambda → stronger shrinkage, more zeros).
%
%        Resampling inference:
%          opts.nPerm         (default 1000) permutations for p-value
%          opts.nBoot         (default 500)  bootstrap resamples for AUC CI
%
%        Learning curve:
%          opts.learningSteps (default 6) number of sample sizes
%
% OUTPUT (results struct)
%   Cross-validated performance (generalization estimate):
%     results.AUC        scalar   mean AUC across repeats×outer folds
%     results.ACC        scalar   mean accuracy across repeats×outer folds
%     results.SENS       scalar   mean sensitivity across repeats×outer folds
%     results.SPEC       scalar   mean specificity across repeats×outer folds
%     results.allAUC     [nRepeats x outerK] fold-level AUC
%     results.allACC     [nRepeats x outerK] fold-level ACC
%     results.allSENS    [nRepeats x outerK] fold-level SENS
%     results.allSPEC    [nRepeats x outerK] fold-level SPEC
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
%     results.AUC_global scalar   AUC of logistic model on mean(X,2)
%
%   Feature stability / interpretability:
%     results.signStability [p x 1] proportion of runs matching mean sign
%     results.selectionFrequency [p x 1] frequency of appearing in topK
%                        absolute weights across runs (topK fixed at 20)
%
%   Permutation test:
%     results.allpermAUC     [nPerm x 1] permuted AUC distribution
%     results.permAUC        scalar mean permuted AUC
%     results.permutation_p  scalar p = mean(permAUC >= observed AUC)
%
%   Bootstrap:
%     results.allbootAUC [nBoot x 1] bootstrap AUC distribution
%     results.bootAUC    scalar mean bootstrap AUC
%     results.AUC_CI     [1 x 2] percentile CI (2.5, 97.5)
%
%   Learning curve:
%     results.learningSizes vector of sample sizes evaluated
%     results.learningAUC   vector of AUC estimates per size
%
% NOTES / INTERPRETATION (high level)
%   - Use results.AUC from nested CV as the primary generalization estimate.
%   - Elastic Net yields sparse solutions (some betas exactly zero), making:
%       featureStability (non-zero proportion) and selectionFrequency
%       particularly informative for interpretation.
%   - Coefficients depend on scaling; this function uses leakage-free z-scoring
%     within each CV fold.
%   - Intercept handling:
%       MATLAB lasso returns an intercept separately (FitInfo.Intercept).
%       This pipeline stores only the p feature coefficients in betaStore/
%       featureWeights. Outer-fold scores are computed here as Xtest*beta
%       (without adding the intercept). This is sufficient for ranking-based
%       metrics such as AUC, but note that the raw scores are not calibrated
%       probabilities.
%
% IMPLEMENTATION NOTES
%   - Inner CV performs a grid search over alphaGrid × lambdaGrid using AUC.
%   - Outer-fold predictions use the full elastic net linear predictor
%       score = X*beta + intercept
%     ensuring consistency with the models evaluated during inner CV.
%   - If some folds lack a class, reduce outerK/innerK.
%
% DEPENDENCIES
%   Requires Statistics and Machine Learning Toolbox:
%     cvpartition, lasso, perfcurve, fitglm, grp2idx
%
% See also: lasso, perfcurve, cvpartition, fitglm

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

rng(1)

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
ACC  = nan(opts.nRepeats,opts.outerK);
SENS = nan(opts.nRepeats,opts.outerK);
SPEC = nan(opts.nRepeats,opts.outerK);

selectedAlpha  = nan(opts.nRepeats,opts.outerK);
selectedLambda = nan(opts.nRepeats,opts.outerK);

interceptStore = nan(opts.nRepeats,opts.outerK);
betaStore = nan(p,opts.outerK,opts.nRepeats);

parfor r = 1:opts.nRepeats

    % local containers for parfor
    AUC_r  = nan(1,opts.outerK);
    ACC_r  = nan(1,opts.outerK);
    SENS_r = nan(1,opts.outerK);
    SPEC_r = nan(1,opts.outerK);

    alpha_r  = nan(1,opts.outerK);
    lambda_r = nan(1,opts.outerK);

    intercept_r = nan(1,opts.outerK);
    beta_r = nan(p,opts.outerK);

    cvOuter = cvpartition(yNum,'KFold',opts.outerK,'Stratify',true);

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

        %% leakage-free scaling

        mu = mean(Xtrain);
        sd = std(Xtrain);
        sd(sd==0)=1;

        Xtrain = (Xtrain-mu)./sd;
        Xtest  = (Xtest-mu)./sd;

        %% inner CV hyperparameter search

        cvInner = cvpartition(ytrain,'KFold',opts.innerK,'Stratify',true);

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

                [B,FitInfo] = lassoglm( ...
                    Xtrain(tr,:),ytr,'binomial',...
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
                            Xtrain,ytrain,'binomial', ...
                            'Alpha',alpha, ...
                            'Lambda',opts.lambdaGrid(l), ...
                            'Standardize',false);

                        bestBeta = Bfull;                    % p x 1
                        bestIntercept = FitInfoFull.Intercept; % scalar

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

        prob = 1./(1+exp(-scores));
        pred = prob > 0.5;

        ACC_r(k) = mean(pred==ytest);

        tp = sum(pred==1 & ytest==1);
        tn = sum(pred==0 & ytest==0);
        fp = sum(pred==1 & ytest==0);
        fn = sum(pred==0 & ytest==1);

        SENS_r(k) = tp/(tp+fn);
        SPEC_r(k) = tn/(tn+fp);

    end

    % write back to shared arrays

    AUC(r,:)  = AUC_r;
    ACC(r,:)  = ACC_r;
    SENS(r,:) = SENS_r;
    SPEC(r,:) = SPEC_r;

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

results.allACC = ACC;
results.ACC = nanmean(ACC(:));

results.allSENS = SENS;
results.SENS = nanmean(SENS(:));

results.allSPEC = SPEC;
results.SPEC = nanmean(SPEC(:));

results.selectedAlpha  = selectedAlpha;
results.selectedLambda = selectedLambda;

results.interceptStore = interceptStore;
results.betaStore = betaStore;
results.featureWeights = featureWeights;

results.featureStability = mean(abs(featureWeights)>0,2);
results.meanFeatureWeight = mean(featureWeights,2);

fprintf('Nested CV AUC = %.3f\n',results.AUC)

%% -------------------------------------------------
% 4. Global signal model
%% -------------------------------------------------

meanBinding = mean(X,2);

mdl = fitglm(meanBinding,yNum,'Distribution','binomial');

scores = predict(mdl,meanBinding);

[~,~,~,AUCg] = perfcurve(yNum,scores,1);

results.AUC_global = AUCg;

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

%% -------------------------------------------------
% 8. Bootstrap CI
%% -------------------------------------------------

bootAUC = nan(opts.nBoot,1);

parfor b=1:opts.nBoot

    idx = randsample(n,n,true);

    bootAUC(b) = quickCV_ENet(X(idx,:),yNum(idx),opts);

end

results.allbootAUC = bootAUC;
results.bootAUC = nanmean(bootAUC);

results.AUC_CI = prctile(bootAUC(~isnan(bootAUC)),[2.5 97.5]);

figure
histogram(bootAUC(~isnan(bootAUC)))
title('Bootstrap AUC')

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