function AUC = bootstrapOOB_PLSDA(X, Y, C, opts)
% bootstrapOOB_PLSDA  Out-of-bag bootstrap AUC for one resample.
%
%   AUC = bootstrapOOB_PLSDA(X, Y, C, opts)
%
%   Draws one bootstrap sample of subjects WITH replacement, trains on it, and
%   evaluates on the out-of-bag subjects. Repeated by the caller nBoot times to
%   build a sampling distribution and a percentile confidence interval.
%
%   Inputs
%     X     [n x p]  features
%     Y     [n x 1]  outcome
%     C     [n x c]  covariates, or []; sliced with the SAME bootstrap index as X
%     opts  struct of pipeline options
%
%   Output
%     AUC   out-of-bag AUC, NaN when the resample was unusable
%
%   Testing out-of-bag rather than on the bootstrap sample itself makes this
%   more conservative than a naive bootstrap, so in small samples it typically
%   sits below the nested-CV estimate. Read it as sampling variability, not as
%   the headline performance number.
%
%   Covariates are resampled with the same index as the features, so the
%   nuisance model is fit on the resampled covariates -- duplicated rows
%   included, since the replicate IS the training set.
%
%   See also plsregress, foldPreprocess, residualizeFold.

n = length(Y);

if numel(unique(Y)) < 2
    AUC = NaN;
    return
end

% Bootstrap sample
idxBoot = randsample(n, n, true);

% OOB = subjects not selected at least once
inBag = false(n,1);
inBag(idxBoot) = true;
oob = ~inBag;

% Need enough OOB cases and both classes present
if sum(oob) < 2
    AUC = NaN;
    return
end

ytrain = Y(idxBoot);
ytest  = Y(oob);

if numel(unique(ytrain)) < 2 || numel(unique(ytest)) < 2
    AUC = NaN;
    return
end

Xtrain = X(idxBoot,:);
Xtest  = X(oob,:);

% leakage-free scaling
Ctrain = []; Ctest = [];
if ~isempty(C)
    Ctrain = C(idxBoot,:);   % duplicated rows: the replicate IS the training set,
    Ctest  = C(oob,:);       % so the nuisance fit uses the resampled covariates
end

[Xtrain, Xtest] = foldPreprocess(Xtrain, Xtest, Ctrain, Ctest, opts.scale);

% cap LV
maxLV = capLV(opts.maxLV, Xtrain);
if maxLV < 1
    AUC = NaN;
    return
end

% inner CV for LV tuning
innerK = min([opts.innerK, floor(length(ytrain)/2), sum(ytrain==0), sum(ytrain==1)]);
if innerK < 2
    AUC = NaN;
    return
end

try
    cvInner = cvpartition(ytrain,'KFold',innerK,'Stratify',true);
catch
    cvInner = cvpartition(length(ytrain),'KFold',innerK);
end

innerAUC = nan(maxLV,1);

for lv = 1:maxLV
    foldAUC = nan(innerK,1);

    for f = 1:innerK
        tr = training(cvInner,f);
        va = test(cvInner,f);

        ytr = ytrain(tr);
        yva = ytrain(va);

        if numel(unique(ytr))<2 || numel(unique(yva))<2
            continue
        end

        [~,~,~,~,beta] = plsregress(Xtrain(tr,:), ytr, lv);
        yhat = [ones(sum(va),1) Xtrain(va,:)] * beta;

        [~,~,~,foldAUC(f)] = perfcurve(yva, yhat, 1);
    end

    innerAUC(lv) = mean(foldAUC,'omitnan');
end

[~,bestLV] = max(innerAUC);

% final fit on bootstrap sample
[~,~,~,~,beta] = plsregress(Xtrain, ytrain, bestLV);

score = [ones(sum(oob),1) Xtest] * beta;
[~,~,~,AUC] = perfcurve(ytest, score, 1);
