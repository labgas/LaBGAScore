function AUC = quickCV_ENet(X,Y,opts)

% quickCV_ENet  Leakage-free quick CV AUC estimate for Elastic Net (ENet).
%
% This helper function provides a *quick* (non-nested) cross-validated AUC
% estimate for Elastic Net models, intended specifically for high-repeat
% diagnostics inside:
%   - permutation testing
%   - bootstrap confidence intervals
%   - learning curves
%
% IMPORTANT (methodological note)
%   This version avoids "test-fold tuning": alpha/lambda are selected using
%   cross-validation *within the training fold only* (via lassoglm with 'CV'),
%   and performance is then evaluated once on the held-out fold.
%
% The function mirrors the main ENet pipelines’ core logic:
%   - leakage-free preprocessing (train-only scaling per fold)
%   - alpha/lambda tuning using training-only CV
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
%        Binary outcome labels. Typically coded as 0/1 (recommended).
%        Folds that end up missing a class are skipped. If Y has <2 unique
%        values overall, the function returns NaN.
%
%   opts struct with fields (required unless stated):
%        opts.outerK      requested outer K for the quick CV estimate.
%                        In this implementation, the actual K is capped to
%                        avoid invalid stratified folds:
%                           K = min([opts.outerK, floor(n/2), minClass])
%                        where minClass is the minority-class count in Y.
%
%        opts.alphaGrid   vector of ENet mixing parameters in [0,1]
%                        alpha=1   lasso; alpha~0 ridge-like.
%                        (May be coarsened internally for speed.)
%
%        opts.lambdaGrid  vector of lambda penalty strengths to evaluate.
%                        (May be coarsened internally for speed.)
%
%        opts.scale       (optional) scaling mode per fold:
%                        'zscore' | 'center' | 'none'
%                        Default is 'zscore' if absent/empty.
%
% OUTPUT
%   AUC  scalar
%        Mean AUC across folds (nanmean over fold AUCs). Returns NaN if:
%        - Y has <2 unique values overall, or
%        - K < 2 after capping, or
%        - all folds are invalid (e.g., missing class after partitioning),
%          or models shrink to all-zero betas in all valid folds.
%
% IMPLEMENTATION NOTES / FIXES VS OLDER VERSIONS
%   1) Training-only hyperparameter tuning:
%      alpha/lambda are chosen using lassoglm with CV on the training fold
%      only, then evaluated on the held-out fold (no information leakage).
%
%   2) Intercept handled consistently:
%      scoring always uses X*B + Intercept (log-odds) for both selection and
%      evaluation.
%
%   3) Robust K selection to prevent stratification warnings:
%      outer K is capped by the minority-class count so each fold can contain
%      both classes when stratifying:
%         K <= minClass.
%      Similarly, inner Kin is capped by the minority-class count *within the
%      training fold*.
%
%   4) Robust cvpartition fallback:
%      if stratified cvpartition fails, falls back to non-stratified KFold.
%
%   5) Speed-oriented coarsening:
%      alphaGrid and lambdaGrid may be downsampled internally for speed; this is
%      appropriate for permutation/bootstrap/learning-curve diagnostics where
%      many repeats are required. For final model selection, use the main
%      nested-CV pipelines.
%
% DEPENDENCIES
%   Requires Statistics and Machine Learning Toolbox:
%     cvpartition, lassoglm, perfcurve
%
% See also: lassoglm, perfcurve, cvpartition, lasso
%
% NOTE
%   This function is designed for speed and diagnostic use. Its AUC is not a
%   replacement for the nested-CV AUC from the main ENet pipeline.

% -----------------------
% Basic checks
% -----------------------
if numel(unique(Y)) < 2
    AUC = NaN; return
end

[n,p] = size(X);

% ---- NEW: cap outer K by minority-class count (avoids stratify warning)
n1 = sum(Y==1);
n0 = sum(Y==0);
minClass = min(n0,n1);

K = min([opts.outerK, floor(n/2), minClass]);
if K < 2
    AUC = NaN; return
end

% robust partition (stratify only if feasible)
try
    cv = cvpartition(Y,'KFold',K,'Stratify',true);
catch
    cv = cvpartition(n,'KFold',K);
end

auc = nan(K,1);

% scaling mode (optional)
scaleMode = 'zscore';
if isfield(opts,'scale') && ~isempty(opts.scale)
    scaleMode = opts.scale;
end

alphaGrid  = opts.alphaGrid;
lambdaGrid = opts.lambdaGrid;

% optional grid coarsening for speed
if numel(alphaGrid) > 4
    alphaGrid = alphaGrid([1 round(end/2) end]);
end
if numel(lambdaGrid) > 12
    lambdaGrid = lambdaGrid(round(linspace(1,numel(lambdaGrid),10)));
end

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

    % leakage-free scaling
    [Xtr, Xte] = applyScaling(Xtr, Xte, scaleMode);

    % ---- NEW: cap inner Kin by minority-class count in TRAIN fold
    n1tr = sum(ytr==1);
    n0tr = sum(ytr==0);
    minClassTr = min(n0tr,n1tr);

    Kin = min([3, floor(numel(ytr)/2), minClassTr]);
    if Kin < 2
        continue
    end

    % robust inner partition
    try
        cvIn = cvpartition(ytr,'KFold',Kin,'Stratify',true);
    catch
        cvIn = cvpartition(length(ytr),'KFold',Kin);
    end

    bestLoss = inf;
    bestB = zeros(p,1);
    bestInt = 0;

    for a = 1:numel(alphaGrid)
        alpha = alphaGrid(a);

        [B,FitInfo] = lassoglm( ...
            Xtr, ytr, 'binomial', ...
            'Alpha', alpha, ...
            'Lambda', lambdaGrid, ...
            'Standardize', false, ...
            'CV', cvIn);

        idx = FitInfo.IndexMinDeviance;
        loss = FitInfo.Deviance(idx);

        if loss < bestLoss
            bestLoss = loss;
            bestB = B(:,idx);
            bestInt = FitInfo.Intercept(idx);
        end
    end

    % if everything shrank to 0, skip fold
    if all(bestB==0)
        continue
    end

    score = Xte*bestB + bestInt; % log-odds
    [~,~,~,auc(k)] = perfcurve(yte, score, 1);
end

AUC = nanmean(auc);

end

% -----------------------
% Local helper: scaling
% -----------------------
function [XtrS, XteS] = applyScaling(Xtr, Xte, mode)
switch lower(mode)
    case 'none'
        XtrS = Xtr; XteS = Xte;
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