function M = quickCV_PLSDA_core(X, Y, C, opts, metric)
% quickCV_PLSDA_core  Fast tuned cross-validated AUC for PLS-DA.
%
%   M = quickCV_PLSDA_core(X, Y, C, opts, metric)
%
%   Single (non-repeated) stratified K-fold cross-validation with train-only
%   preprocessing and inner-CV selection of the number of latent variables.
%   Shared implementation behind quickCV_PLSDA (ROC AUC) and quickCV_PLSDA_PR
%   (precision-recall AUC), which previously existed as two files differing by
%   two lines.
%
%   Inputs
%     X       [n x p]  features
%     Y       [n x 1]  binary 0/1 labels
%     C       [n x c]  covariates sliced with the SAME index as X, or []
%     opts    struct; reads .outerK, .innerK, .maxLV, .scale
%     metric  'roc' (default) | 'pr'
%
%   Output
%     M       mean AUC across folds, NaN if no fold was usable
%
%   Guards this adds over the previous implementation
%     - K is capped by the minority-class count as well as floor(n/2), matching
%       quickCV_ENet. Without it, a stratified partition asking for more folds
%       than there are minority-class members emits a warning and silently
%       produces folds missing a class, so every fold is skipped and the result
%       is NaN.
%     - An explicit K < 2 return. Previously cvpartition(...,'KFold',1) threw,
%       the catch block called cvpartition again with the same bad K, and the
%       error escaped and aborted the whole enclosing parfor.
%
%   Latent variables are tuned here so this estimator matches the one producing
%   the observed statistic; see quickCV_PLSR for why that matters.
%
%   See also quickCV_PLSDA, quickCV_PLSDA_PR, foldPreprocess, capLV.

if nargin < 5 || isempty(metric)
    metric = 'roc';
end

M = NaN;

if numel(unique(Y)) < 2
    return
end

n = numel(Y);

n1 = sum(Y == max(Y));
n0 = n - n1;
minClass = min(n0, n1);

K = min([opts.outerK, floor(n/2), minClass]);

if K < 2
    return
end

if ~isempty(C) && size(C,1) ~= n
    error('quickCV_PLSDA_core:covariateRows', ...
        'C has %d rows but Y has %d. Slice covariates with the same index as X.', size(C,1), n);
end

try
    cv = cvpartition(Y,'KFold',K,'Stratify',true);
catch
    cv = cvpartition(n,'KFold',K);
end

vals = nan(K,1);

for k = 1:K

    tr = training(cv,k);
    te = test(cv,k);

    ytr = Y(tr);
    yte = Y(te);

    if numel(unique(ytr)) < 2 || numel(unique(yte)) < 2
        continue
    end

    Ctr = []; Cte = [];
    if ~isempty(C), Ctr = C(tr,:); Cte = C(te,:); end

    [Xtr, Xte] = foldPreprocess(X(tr,:), X(te,:), Ctr, Cte, opts.scale);

    maxLV = capLV(opts.maxLV, Xtr);
    if maxLV < 1
        continue
    end

    % ---- inner CV selection of LV, training fold only ----
    innerK = min([opts.innerK, numel(ytr)-1, min(sum(ytr==1), sum(ytr==0))]);
    if innerK < 2
        continue
    end

    try
        cvIn = cvpartition(ytr,'KFold',innerK,'Stratify',true);
    catch
        cvIn = cvpartition(numel(ytr),'KFold',innerK);
    end

    Xtr_raw = X(tr,:);
    innerVal = nan(maxLV,1);

    for lv = 1:maxLV
        foldVal = nan(innerK,1);
        for f = 1:innerK
            t2 = training(cvIn,f);
            v2 = test(cvIn,f);
            y2 = ytr(t2);
            yv = ytr(v2);
            if numel(unique(y2)) < 2 || numel(unique(yv)) < 2
                continue
            end
            C2 = []; Cv = [];
            if ~isempty(Ctr), C2 = Ctr(t2,:); Cv = Ctr(v2,:); end
            [X2, Xv] = foldPreprocess(Xtr_raw(t2,:), Xtr_raw(v2,:), C2, Cv, opts.scale);
            if capLV(lv, X2) < lv
                continue
            end
            [~,~,~,~,b] = plsregress(X2, y2, lv);
            sc = [ones(sum(v2),1) Xv] * b;
            foldVal(f) = scoreMetric(yv, sc, metric);
        end
        innerVal(lv) = mean(foldVal,'omitnan');
    end

    if all(isnan(innerVal))
        continue
    end
    [~,bestLV] = max(innerVal);

    [~,~,~,~,beta] = plsregress(Xtr, ytr, bestLV);
    score = [ones(sum(te),1) Xte] * beta;

    vals(k) = scoreMetric(yte, score, metric);

end

M = mean(vals,'omitnan');

end

% -------------------------------------------------------------------------
function v = scoreMetric(ytrue, score, metric)
v = NaN;
if numel(unique(score)) < 2 || ~all(isfinite(score))
    return
end
try
    if strcmpi(metric,'pr')
        [~,~,~,v] = perfcurve(ytrue, score, 1, 'xCrit','reca', 'yCrit','prec');
    else
        [~,~,~,v] = perfcurve(ytrue, score, 1);
    end
catch
    v = NaN;
end
end
