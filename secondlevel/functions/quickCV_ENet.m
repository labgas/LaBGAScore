function AUC = quickCV_ENet(X,Y,opts)

% Leakage-free quick CV AUC estimate for Elastic Net (ENet).
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
% In addition, this version includes robustness safeguards designed for
% large permutation/bootstrap runs in small-sample settings:
%   - outer K is capped by the minority-class count to avoid invalid folds
%   - inner Kin is capped similarly within the training fold
%   - folds with missing class or invalid AUC estimates are skipped (NaN)
%     rather than forced to 0.5 to prevent artificial spikes in permutation
%     distributions
%   - intercept-only fallback is used if no usable model is found
%   - diagnostic counters summarizing fold outcomes are printed
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
%        Folds that end up missing a class are skipped (NaN). If Y has <2
%        unique values overall, the function returns NaN.
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
%        opts.scale       (optional) scaling mode per fold:
%                        'zscore' | 'center' | 'none'
%                        Default is 'zscore' if absent/empty.
%
% OUTPUT
%   AUC  scalar
%        Mean AUC across folds (nanmean over fold AUCs).
%        Returns NaN if:
%        - Y has <2 unique values overall, or
%        - K < 2 after capping.
%
%        Returns 0.5 if all folds are invalid (all fold AUC values are NaN).
%
% IMPLEMENTATION NOTES / FIXES VS OLDER VERSIONS
%   1) Training-only hyperparameter tuning:
%      alpha/lambda are chosen using lassoglm with CV on the training fold
%      only, then evaluated on the held-out fold (no information leakage).
%
%   2) Intercept handled consistently:
%      scoring always uses X*B + Intercept (log-odds). Intercept-only fallback
%      uses logit(mean(ytr)).
%
%   3) Robust K selection to prevent stratification warnings:
%      outer K is capped by the minority-class count so each fold can contain
%      both classes when stratifying:
%         K <= minClass.
%      Similarly, inner Kin is capped by the minority-class count within the
%      training fold.
%
%   4) Robust lambda selection:
%      lambda index preference is Index1SE (stability). If Index1SE is invalid
%      or produces an all-zero model, the algorithm falls back to
%      IndexMinDeviance. If that still produces all-zero coefficients, the
%      best non-zero solution along the lambda path (among finite deviance
%      entries) is selected if available.
%
%   5) Hardcoded lassoglm parameters (for stability in small-n settings):
%         'NumLambda'   = 25
%         'LambdaRatio' = 1e-3
%         'MaxIter'     = 1e4
%         'RelTol'      = 1e-3
%
%      These settings limit extremely small lambda values that often cause
%      convergence problems in high p >> n scenarios while keeping the
%      solution path sufficiently rich for permutation diagnostics.
%
%   6) Diagnostic counters:
%      A diagnostic struct (diag) is printed summarizing fold outcomes:
%         nFolds
%         nFoldMissingClass
%         nKinTooSmall
%         nNoFitAllAlpha
%         nAllZeroChosen
%         nPerfcurveFail
%         nPerfcurveNonFinite
%         nValidFolds
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

diag = struct;
diag.nFolds = K;
diag.nFoldMissingClass = 0;
diag.nKinTooSmall = 0;
diag.nNoFitAllAlpha = 0;
diag.nAllZeroChosen = 0;
diag.nPerfcurveFail = 0;
diag.nPerfcurveNonFinite = 0;
diag.nValidFolds = 0;

% scaling mode (optional)
scaleMode = 'zscore';
if isfield(opts,'scale') && ~isempty(opts.scale)
    scaleMode = opts.scale;
end

alphaGrid  = opts.alphaGrid;

% optional grid coarsening for speed
if numel(alphaGrid) > 4
    alphaGrid = alphaGrid([1 round(end/2) end]);
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
        diag.nFoldMissingClass = diag.nFoldMissingClass + 1;
        auc(k) = NaN; % skip invalid fold (prevents 0.5 pile-up)
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
        diag.nKinTooSmall = diag.nKinTooSmall + 1;
        % fall back to intercept-only evaluation on this fold
        bestInt = logitSafe(mean(ytr));
        score = bestInt * ones(sum(te),1);
        [~,~,~,auc(k)] = perfcurve(yte, score, 1);
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
    bestInt = logitSafe(mean(ytr));

    anyFit = false;
    
    for a = 1:numel(alphaGrid)
        alpha = alphaGrid(a);
        
        % Train-only tuning for lambda using CV on training data
        % Use the SAME lambdaGrid we want to consider.
        try
        [B,FitInfo] = lassoglm( ...
            Xtr, ytr, 'binomial', ...
            'Alpha', alpha, ...
            'Standardize', false, ...
            'CV', cvIn,...
            'NumLambda', 25, ...
            'LambdaRatio', 1e-3, ...   % try 1e-2 or 1e-3; avoid 1e-4+ in tiny n
            'MaxIter', 1e4,...
            'RelTol', 1e-3);
        anyFit = true;
        catch
            continue
        end

        % Choose lambda index (start with 1SE for stability)
        idx1 = FitInfo.Index1SE;
        idxm = FitInfo.IndexMinDeviance;

        idx = idx1;

        % If 1SE is invalid or all-zero, try min deviance
        if isempty(idx) || idx<1 || idx>size(B,2) || all(B(:,idx)==0)
            idx = idxm;
        end

        % If still all-zero, choose the best NON-zero model on the path (if any)
        if ~isempty(idx) && idx>=1 && idx<=size(B,2) && all(B(:,idx)==0)
            nonzeroCols = find(any(B~=0,1));        % lambdas with at least 1 active feature
            finiteDev   = find(isfinite(FitInfo.Deviance(:))');
            candidates  = intersect(nonzeroCols, finiteDev);
        
            if ~isempty(candidates)
                [~,ii] = min(FitInfo.Deviance(candidates));
                idx = candidates(ii);
            end
        end

        % Robustify index one last time
        dev = FitInfo.Deviance;
        if isempty(idx) || isnan(idx) || idx<1 || idx>numel(dev) || ~isfinite(dev(idx))
            finiteIdx = find(isfinite(dev));
            if isempty(finiteIdx)
                continue
            end
            [~,m] = min(dev(finiteIdx));
            idx = finiteIdx(m);
        end

        loss = dev(idx);

        if loss < bestLoss
            bestLoss = loss;
            bestB = B(:,idx);
            bestInt = FitInfo.Intercept(idx);
        end
    end

    % IMPORTANT CHANGE:
    % Do NOT skip all-zero solutions; evaluate intercept-only model instead.
    if ~anyFit
        diag.nNoFitAllAlpha = diag.nNoFitAllAlpha + 1;
        % no alpha worked -> intercept-only
        bestInt = logitSafe(mean(ytr));
        score = bestInt * ones(sum(te),1);
    else
        if all(bestB==0)
            diag.nAllZeroChosen = diag.nAllZeroChosen + 1;
            bestInt = logitSafe(mean(ytr));
            score = bestInt * ones(sum(te),1);
        else
            score = Xte*bestB + bestInt; % if bestB==0 => constant score
        end
    end

    try
    [~,~,~,auc(k)] = perfcurve(yte, score, 1);
        if ~isfinite(auc(k))
            diag.nPerfcurveNonFinite = diag.nPerfcurveNonFinite + 1;
            auc(k) = NaN; 
        end
    catch
        diag.nPerfcurveFail = diag.nPerfcurveFail + 1;
        auc(k) = NaN;
    end

end

diag.nValidFolds = sum(isfinite(auc));
disp(diag);

if all(isnan(auc))
    AUC = 0.5;
else
    AUC = nanmean(auc);
end

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

function z = logitSafe(p)
% stable logit for p in [0,1]
p = min(max(p, 1e-6), 1-1e-6);
z = log(p/(1-p));
end