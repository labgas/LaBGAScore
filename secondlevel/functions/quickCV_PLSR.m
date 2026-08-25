function Q2 = quickCV_PLSR(X, Y, C, opts)
% quickCV_PLSR  Fast tuned cross-validated Q2 for PLSR, used for resampling.
%
%   Q2 = quickCV_PLSR(X, Y, C, opts)
%
%   Single (non-repeated) K-fold cross-validation with train-only preprocessing
%   and inner-CV selection of the number of latent variables. Used by the
%   permutation, learning-curve and matched-observed stages of
%   PLSR_neuroimaging_pipeline, where the full repeated nested CV would be too
%   expensive to run thousands of times.
%
%   Inputs
%     X     [n x p]  features
%     Y     [n x 1]  continuous outcome
%     C     [n x c]  covariates sliced with the SAME index as X, or []
%     opts  struct; reads .outerK, .innerK, .maxLV, .scale
%
%   Output
%     Q2    mean out-of-sample Q2 across folds, NaN if no fold was usable
%
%   Why this selects latent variables
%     This function supplies the permutation null, while the observed statistic
%     comes from repeated nested CV WITH inner tuning. When the null used the
%     maximum LV count and the observed statistic used a tuned one, the two
%     sides of the p-value were different estimators: the untuned estimator
%     overfits permuted data far harder, pushing the null well below the
%     observed value. Measured on true-null data, that mismatch produced a
%     false-positive rate near 40% at alpha = 0.05. Tuning here, and comparing
%     against a matched observed statistic computed by this same function,
%     restores calibration.
%
%   See also PLSR_neuroimaging_pipeline, foldPreprocess, capLV.

if nargin < 4
    error('quickCV_PLSR:nargin', ...
        'quickCV_PLSR now requires (X, Y, C, opts). Pass C = [] when no covariates are used.');
end

n = numel(Y);
K = min(opts.outerK, floor(n/2));

if K < 2
    Q2 = NaN;
    return
end

if ~isempty(C) && size(C,1) ~= n
    error('quickCV_PLSR:covariateRows', ...
        'C has %d rows but Y has %d. Slice covariates with the same index as X.', size(C,1), n);
end

cv = cvpartition(n,'KFold',K);
q2 = nan(K,1);

for k = 1:K

    tr = training(cv,k);
    te = test(cv,k);

    ytr = Y(tr);
    yte = Y(te);

    if numel(ytr) < 3 || numel(yte) < 2
        continue
    end

    Ctr = []; Cte = [];
    if ~isempty(C), Ctr = C(tr,:); Cte = C(te,:); end

    [Xtr, Xte] = foldPreprocess(X(tr,:), X(te,:), Ctr, Cte, opts.scale);

    maxLV = capLV(opts.maxLV, Xtr);
    if maxLV < 1
        continue
    end

    % ---- inner CV selection of LV, on the training fold only ----
    innerK = min([opts.innerK, numel(ytr)-1]);
    if innerK < 2
        continue
    end

    cvIn = cvpartition(numel(ytr),'KFold',innerK);
    innerQ2 = nan(maxLV,1);

    Xtr_raw = X(tr,:);

    for lv = 1:maxLV
        foldQ2 = nan(innerK,1);
        for f = 1:innerK
            t2 = training(cvIn,f);
            v2 = test(cvIn,f);
            y2 = ytr(t2);
            yv = ytr(v2);
            if numel(y2) < 3 || numel(yv) < 2
                continue
            end
            C2 = []; Cv = [];
            if ~isempty(Ctr), C2 = Ctr(t2,:); Cv = Ctr(v2,:); end
            [X2, Xv] = foldPreprocess(Xtr_raw(t2,:), Xtr_raw(v2,:), C2, Cv, opts.scale);
            if capLV(lv, X2) < lv
                continue
            end
            [~,~,~,~,b] = plsregress(X2, y2, lv);
            yh = [ones(sum(v2),1) Xv] * b;
            d  = sum((yv - mean(y2)).^2);
            if d <= 0
                continue
            end
            foldQ2(f) = 1 - sum((yv - yh).^2) / d;
        end
        innerQ2(lv) = mean(foldQ2,'omitnan');
    end

    if all(isnan(innerQ2))
        continue
    end
    [~,bestLV] = max(innerQ2);

    [~,~,~,~,beta] = plsregress(Xtr, ytr, bestLV);
    yhat = [ones(sum(te),1) Xte] * beta;

    denom = sum((yte - mean(ytr)).^2);
    if denom <= 0
        continue
    end

    q2(k) = 1 - sum((yte - yhat).^2) / denom;

end

Q2 = mean(q2,'omitnan');

end
