function results = TSPO_PLSDA_pipeline(X,Y,opts)

% TSPO_PLSDA_pipeline  Robust PLS-DA pipeline for neuroimaging ROI features.
%
% This function implements Partial Least Squares Discriminant Analysis (PLS-DA)
% for binary classification with repeated nested k-fold cross-validation.
% It is designed for neuroimaging feature matrices (subjects × features) such
% as PET ROI binding, fMRI ROI betas, morphometry, or connectivity-derived
% measures. The architecture emphasizes:
%   - leakage-free preprocessing (train-only scaling in every fold)
%   - inner CV selection of number of latent variables (LVs)
%   - outer CV estimation of generalization performance
%   - resampling-based inference (permutation p-value, bootstrap CI)
%   - stability/interpretability metrics (VIP, stabilityZ, sign stability,
%     top-K selection frequency)
%   - learning curve (AUC vs sample size)
%
% USAGE
%   results = TSPO_PLSDA_pipeline(X, Y)
%   results = TSPO_PLSDA_pipeline(X, Y, opts)
%
% INPUTS
%   X    [n x p] numeric
%        Feature matrix (n subjects, p features). Each column is a feature
%        (e.g., ROI binding, ROI beta, graph metric). Missing values should
%        be handled upstream (impute/remove) prior to calling this function.
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
%     results.selectedLV [nRepeats x outerK] selected LV per outer fold
%     results.betaStore  [p+1 x outerK x nRepeats] PLS regression betas
%                        (row 1 is intercept; rows 2..p+1 correspond to features)
%     results.featureWeights [p x (nRepeats*outerK)] betas (no intercept),
%                        stacked across all outer folds/repeats
%     results.meanFeatureWeight [p x 1] mean featureWeights across runs
%     results.featureStability  [p x 1] proportion of runs where |beta|>0
%                        (in PLS-DA this is typically ~1 because betas are
%                         rarely exactly zero; included for plot symmetry)
%
%   Global baseline (interpretation only):
%     results.AUC_global scalar   AUC of logistic model on mean(X,2)
%
%   Final model on all data (interpretation only; NOT for performance):
%     results.finalLV        scalar median selected LV across all runs
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
%   - Bootstrap AUC can be higher than nested CV AUC (often optimistic,
%     especially for small n) and is best viewed as sampling variability.
%   - VIP and stabilityZ provide complementary interpretability:
%       VIP: importance in the final fitted model
%       stabilityZ: robustness of the effect across CV runs
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
% 2. Repeated Nested Cross-Validation
%% -------------------------------------------------

AUC = nan(opts.nRepeats,opts.outerK);
ACC = nan(opts.nRepeats,opts.outerK);
SENS = nan(opts.nRepeats,opts.outerK);
SPEC = nan(opts.nRepeats,opts.outerK);

selectedLV = nan(opts.nRepeats,opts.outerK);
betaStore = nan(p+1,opts.outerK,opts.nRepeats);
featureWeights = nan(p, opts.nRepeats*opts.outerK);

for r = 1:opts.nRepeats

    cvOuter = cvpartition(yNum,'KFold',opts.outerK,'Stratify',true);

    for k = 1:opts.outerK

        trainIdx = training(cvOuter,k);
        testIdx = test(cvOuter,k);

        ytrain = yNum(trainIdx);
        ytest = yNum(testIdx);

        % skip if class missing
        if numel(unique(ytrain))<2 || numel(unique(ytest))<2
            continue
        end

        Xtrain = X(trainIdx,:);
        Xtest = X(testIdx,:);

        %% leakage-free scaling

        mu = mean(Xtrain);
        sd = std(Xtrain);
        sd(sd==0)=1;

        Xtrain = (Xtrain-mu)./sd;
        Xtest  = (Xtest-mu)./sd;

        %% inner CV LV tuning

        maxLV = min([opts.maxLV rank(Xtrain)-1 size(Xtrain,1)-2]);

        innerAUC = nan(maxLV,1);
        cvInner = cvpartition(ytrain,'KFold',opts.innerK,'Stratify',true);

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

                [~,~,~,~,beta] = plsregress(Xtrain(tr,:),ytr,lv);

                yhat = [ones(sum(va),1) Xtrain(va,:)]*beta;

                [~,~,~,foldAUC(f)] = perfcurve(yva,yhat,1);

            end

            innerAUC(lv) = nanmean(foldAUC);

        end

        [~,bestLV] = max(innerAUC);
        selectedLV(r,k) = bestLV;

        %% fit model

        [~,~,~,~,beta] = plsregress(Xtrain,ytrain,bestLV);
        betaStore(:,k,r) = beta;
        featureWeights(:,(r-1)*opts.outerK+k) = beta(2:end);

        scores = [ones(sum(testIdx),1) Xtest]*beta;

        [~,~,~,AUC(r,k)] = perfcurve(ytest,scores,1);

        pred = scores>0.5;

        ACC(r,k) = mean(pred==ytest);

        tp = sum(pred==1 & ytest==1);
        tn = sum(pred==0 & ytest==0);
        fp = sum(pred==1 & ytest==0);
        fn = sum(pred==0 & ytest==1);

        SENS(r,k) = tp/(tp+fn);
        SPEC(r,k) = tn/(tn+fp);

    end
end

results.allAUC = AUC;
results.AUC  = nanmean(AUC(:));
results.allACC = ACC;
results.ACC  = nanmean(ACC(:));
results.allSENS = SENS;
results.SENS = nanmean(SENS(:));
results.allSPEC = SPEC;
results.SPEC = nanmean(SPEC(:));

results.selectedLV = selectedLV;
results.betaStore = betaStore;
results.featureWeights = featureWeights;

results.featureStability = mean(abs(featureWeights) > 0,2);
results.meanFeatureWeight = mean(featureWeights,2);

fprintf('Nested CV AUC = %.3f\n',results.AUC)

%% -------------------------------------------------
% 3. Global signal model
%% -------------------------------------------------

meanBinding = mean(X,2);

mdl = fitglm(meanBinding,yNum,'Distribution','binomial');
scores = predict(mdl,meanBinding);

[~,~,~,AUCg] = perfcurve(yNum,scores,1);

results.AUC_global = AUCg;

%% -------------------------------------------------
% 4. Final model (interpretation only)
%% -------------------------------------------------

mu = mean(X);
sd = std(X);
sd(sd==0)=1;

Xz = (X-mu)./sd;

finalLV = round(nanmedian(selectedLV(:)));

[XL,YL,XS,YS,beta,PCTVAR,MSE,stats] = plsregress(Xz,yNum,finalLV);

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
meanBeta = nanmean(betaMat,2);
sdBeta = nanstd(betaMat,[],2);

results.meanBeta = meanBeta;
results.sdBeta = sdBeta;
results.stabilityZ = meanBeta ./ sdBeta;

%% -------------------------------------------------
% 7. Permutation testing
%% -------------------------------------------------

permAUC = nan(opts.nPerm,1);

parfor i=1:opts.nPerm
    yp = yNum(randperm(n));
    permAUC(i) = quickCV(X,yp,opts);
end

results.allpermAUC = permAUC;
results.permAUC = nanmean(permAUC);
results.permutation_p = mean(permAUC >= results.AUC,'omitnan');

figure
histogram(permAUC(~isnan(permAUC)))
hold on
xline(results.AUC)
title('Permutation AUC distribution')

%% -------------------------------------------------
% 8. Bootstrap AUC CI
%% -------------------------------------------------

bootAUC = nan(opts.nBoot,1);

parfor b=1:opts.nBoot
    idx = randsample(n,n,true);
    bootAUC(b) = quickCV(X(idx,:),yNum(idx),opts);
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

    lcAUC(i) = quickCV(X(idx,:),yNum(idx),opts);

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

signStability = mean(sign(results.meanFeatureWeight) == ...
    sign(results.betaStore(2:end,:,:,:)),[3 4]);
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

for i = 1:size(results.featureWeights,2)
    [~,idx] = sort(abs(results.featureWeights(:,i)),'descend');
    freq(idx(1:topK)) = freq(idx(1:topK)) + 1;
end

freq = freq / size(results.featureWeights,2);

results.selectionFrequency = freq;


end

%% inline functions

function AUC = quickCV(X,Y,opts)

if numel(unique(Y)) < 2
    AUC = NaN;
    return
end

n = length(Y);
K = min(opts.outerK,floor(n/2));

cv = cvpartition(Y,'KFold',K,'Stratify',true);

auc = nan(K,1);

for k=1:K

    tr = training(cv,k);
    te = test(cv,k);

    ytr = Y(tr);
    yte = Y(te);

    if numel(unique(ytr))<2 || numel(unique(yte))<2
        continue
    end

    Xtr = X(tr,:);
    Xte = X(te,:);

    mu = mean(Xtr);
    sd = std(Xtr);
    sd(sd==0)=1;

    Xtr = (Xtr-mu)./sd;
    Xte = (Xte-mu)./sd;

    lv = min([opts.maxLV rank(Xtr)-1 size(Xtr,1)-2]);

    [~,~,~,~,beta] = plsregress(Xtr,ytr,lv);

    score = [ones(sum(te),1) Xte]*beta;

    [~,~,~,auc(k)] = perfcurve(yte,score,1);

end

AUC = nanmean(auc);

end