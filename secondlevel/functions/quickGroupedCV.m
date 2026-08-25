function AUC = quickGroupedCV(X, Y, subjIdx, C, K, opts)
% quickGroupedCV  Subject-blocked grouped CV for the paired PLS-DA permutation null.
%
%   AUC = quickGroupedCV(X, Y, subjIdx, C, K, opts)
%
%   Folds are formed over SUBJECTS, so both observations of a subject always
%   land on the same side of the split. Pass C = [] when no covariates are used;
%   C must be sliced with the same row index as X.
%
%   See also PLSDA_paired_neuroimaging_pipeline, makeGroupedFolds, foldPreprocess.

nSubj = max(subjIdx);
K = min(K, nSubj);
foldID = makeGroupedFolds(nSubj, K);

auc = nan(K,1);

for k = 1:K
    teSubj = foldID == k;
    trSubj = ~teSubj;

    tr = trSubj(subjIdx);
    te = teSubj(subjIdx);

    ytr = Y(tr);
    yte = Y(te);

    Ctr = []; Cte = [];
    if ~isempty(C), Ctr = C(tr,:); Cte = C(te,:); end

    if numel(unique(ytr)) < 2 || numel(unique(yte)) < 2
        continue
    end

    Xtr = X(tr,:);
    Xte = X(te,:);
    [Xtr, Xte] = foldPreprocess(Xtr, Xte, Ctr, Cte, opts.scale);

    lv = capLV(opts.maxLV, Xtr);
    if lv < 1
        continue
    end

    [~,~,~,~,beta] = plsregress(Xtr, ytr, lv);
    score = [ones(sum(te),1) Xte] * beta;
    [~,~,~,auc(k)] = perfcurve(yte, score, 1);
end

AUC = mean(auc,'omitnan');
end
