function AUC = bootstrapOOB_ENet(X, Y, C, opts)
% bootstrapOOB_ENet  Out-of-bag bootstrap AUC for one resample.
%
%   AUC = bootstrapOOB_ENet(X, Y, C, opts)
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
%   See also lassoglm, foldPreprocess, residualizeFold.

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
Ctrain = []; Ctest = [];
if ~isempty(C)
    Ctrain = C(idxBoot,:);   % duplicated rows: the replicate IS the training set,
    Ctest  = C(oob,:);       % so the nuisance fit uses the resampled covariates
end

[Xtrain, Xtest] = foldPreprocess(Xtrain, Xtest, Ctrain, Ctest, opts.scale);

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
            'Lambda', enetLambdaGrid(opts), ...
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

    meanAUC = mean(foldAUC,1,'omitnan');

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
    seAUC = std(foldAUC,[],1,'omitnan') ./ sqrt(sum(~isnan(foldAUC),1));
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

