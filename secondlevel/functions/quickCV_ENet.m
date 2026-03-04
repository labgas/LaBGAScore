function AUC = quickCV_ENet(X,Y,opts)

% quickCV_ENet  Fast cross-validated AUC estimate for Elastic Net (ENet).
%
% This helper function provides a *quick* (non-nested) cross-validated AUC
% estimate for Elastic Net models, and is intended for use inside:
%   - permutation testing
%   - bootstrap confidence intervals
%   - learning curves
%
% The function mirrors the main ENet pipelines’ core logic:
%   - leakage-free preprocessing (train-only scaling per fold)
%   - alpha/lambda grid search within each CV fold (coarsened for speed)
%   - scoring uses the stored intercept consistently (X*B + Intercept)
%
% USAGE
%   AUC = quickCV_ENet(X, Y, opts)
%
% INPUTS
%   X    [n x p] numeric
%        Feature matrix (n subjects, p features). Each column is a feature
%        (e.g., ROI binding, ROI beta, edge weight, graph metric).
%        Missing values should be handled upstream (impute/remove).
%
%   Y    [n x 1] labels
%        Binary outcome labels (0/1, logical, categorical, string, etc.).
%        Must contain two classes in the data used for a fold; otherwise that
%        fold is skipped. If Y has <2 unique values overall, returns NaN.
%
%   opts struct with fields (required unless stated):
%        opts.outerK      outer K for the quick CV estimate
%                        (K is set to min(opts.outerK, floor(n/2)))
%
%        opts.alphaGrid   vector of ENet mixing parameters in [0,1]
%                        alpha=1   lasso; alpha~0 ridge-like
%
%        opts.lambdaGrid  vector of lambda penalty strengths to evaluate
%
%        opts.scale       (optional) scaling mode per fold:
%                        'zscore' | 'center' | 'none'
%                        Default is 'zscore' if absent/empty.
%
% OUTPUT
%   AUC  scalar
%        Mean AUC across folds (nanmean over fold AUCs). Returns NaN if:
%        - Y has <2 unique values overall, or
%        - K < 2, or
%        - all folds are invalid (e.g., missing class after partitioning),
%          or models shrink to all-zero betas in all valid folds.
%
% IMPLEMENTATION NOTES / FIXES VS OLDER VERSIONS
%   1) Lambda indexing alignment:
%      lasso is fit with the SAME (possibly downsampled) lambdaGrid that is
%      evaluated, preventing misalignment between B(:,l) and lambda.
%
%   2) Intercept handled consistently:
%      scoring always uses X*B + Intercept (both for selection and final AUC).
%
%   3) Correct all-zero handling:
%      fold is skipped if bestBeta is all zeros (instead of checking a stale
%      loop variable).
%
%   4) Robust cvpartition fallback:
%      if stratified cvpartition fails, falls back to non-stratified KFold.
%
%   5) Speed-oriented coarsening:
%      alphaGrid and lambdaGrid are downsampled internally for speed; this is
%      appropriate for permutation/bootstrap/learning-curve diagnostics where
%      many repeats are required. For final model selection, use the main
%      nested-CV pipelines.
%
% DEPENDENCIES
%   Requires Statistics and Machine Learning Toolbox:
%     cvpartition, lasso, perfcurve
%
% See also: lasso, perfcurve, cvpartition
%
% NOTE
%   This function is designed for speed and diagnostic use. Its AUC is not a
%   replacement for the nested-CV AUC from the main pipeline.

% -----------------------
% Pull grids and coarsen
% -----------------------
alphaGrid  = opts.alphaGrid;
lambdaGrid = opts.lambdaGrid;

% speedups: coarsen grids (keep endpoints + middle for alpha; ~10 points for lambda)
if numel(alphaGrid) > 4
    alphaGrid = alphaGrid([1 round(end/2) end]);
end

if numel(lambdaGrid) > 12
    lambdaGrid = lambdaGrid(round(linspace(1,numel(lambdaGrid),10)));
end

% -----------------------
% Basic checks
% -----------------------
if numel(unique(Y)) < 2
    AUC = NaN;
    return
end

[n,p] = size(X);

K = min(opts.outerK, floor(n/2));
if K < 2
    AUC = NaN;
    return
end

% robust partition
try
    cv = cvpartition(Y,'KFold',K,'Stratify',true);
catch
    cv = cvpartition(n,'KFold',K);
end

auc = nan(K,1);

% -----------------------
% CV loop
% -----------------------
for k = 1:K

    tr = training(cv,k);
    te = test(cv,k);

    ytr = Y(tr);
    yte = Y(te);

    if numel(unique(ytr))<2 || numel(unique(yte))<2
        continue
    end

    Xtr = X(tr,:);
    Xte = X(te,:);

    % leakage-free scaling (generic)
    scaleMode = 'zscore';
    if isfield(opts,'scale') && ~isempty(opts.scale)
        scaleMode = opts.scale;
    end
    [Xtr, Xte] = applyScaling(Xtr, Xte, scaleMode);

    bestAUC = -inf;
    bestBeta = zeros(p,1);
    bestIntercept = 0;

    for a = 1:numel(alphaGrid)

        alpha = alphaGrid(a);

        % IMPORTANT: fit lasso with the SAME lambdaGrid we evaluate
        [B,FitInfo] = lasso( ...
            Xtr,ytr, ...
            'Alpha',alpha, ...
            'Lambda',lambdaGrid, ...
            'Standardize',false);

        for l = 1:numel(lambdaGrid)

            score = Xte*B(:,l) + FitInfo.Intercept(l);
            [~,~,~,aucTemp] = perfcurve(yte,score,1);

            if aucTemp > bestAUC
                bestAUC = aucTemp;
                bestBeta = B(:,l);
                bestIntercept = FitInfo.Intercept(l);
            end

        end

    end

    % if model shrank everything to 0, skip fold
    if all(bestBeta==0)
        continue
    end

    score = Xte*bestBeta + bestIntercept;
    [~,~,~,auc(k)] = perfcurve(yte,score,1);

end

AUC = nanmean(auc);

end

% -----------------------
% Local helper: scaling
% -----------------------
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