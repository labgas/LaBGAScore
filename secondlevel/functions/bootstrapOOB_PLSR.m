function Q2 = bootstrapOOB_PLSR(X, Y, C, opts)
% bootstrapOOB_PLSR  Out-of-bag bootstrap Q2 for one resample.
%
%   Q2 = bootstrapOOB_PLSR(X, Y, C, opts)
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
%     Q2   out-of-bag Q2, NaN when the resample was unusable
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
Ctrain = []; Ctest = [];
if ~isempty(C)
    Ctrain = C(idxBoot,:);   % duplicated rows: the replicate IS the training set,
    Ctest  = C(oob,:);       % so the nuisance fit uses the resampled covariates
end

[Xtrain, Xtest] = foldPreprocess(Xtrain, Xtest, Ctrain, Ctest, opts.scale);

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

    innerQ2(lv) = mean(foldQ2,'omitnan');
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

