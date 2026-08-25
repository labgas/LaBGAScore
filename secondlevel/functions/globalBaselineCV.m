function out = globalBaselineCV(X, Y, C, opts, kind)
% globalBaselineCV  Cross-validated global-signal baseline.
%
%   out = globalBaselineCV(X, Y, C, opts, kind)
%
%   Companion to the in-sample global-signal baseline that every pipeline
%   already reports. The in-sample version fits fitlm/fitglm on all subjects
%   and predicts those same subjects, so results.Q2_global / results.AUC_global
%   are IN-SAMPLE numbers sitting next to a CROSS-VALIDATED results.Q2 /
%   results.AUC. Comparing them directly flatters the baseline.
%
%   This function reruns that baseline through the same repeated outer CV the
%   model uses, so the two are finally comparable. Both are reported: the
%   in-sample fields keep their previous meaning and values, and the *_cv
%   fields added alongside them are the ones to compare against the model.
%
%   Inputs
%     X     [n x p]  features
%     Y     [n x 1]  outcome (continuous for 'regress', binary 0/1 otherwise)
%     C     [n x c]  covariates, or []
%     opts  struct; reads .outerK, .nRepeats, .globalFun
%     kind  'regress' | 'classify'
%
%   Output (struct)
%     'regress'  : .Q2
%     'classify' : .AUC, .AUC_PR
%
%   The global feature is a per-subject row summary (mean or median across
%   features, or opts.globalFun), so computing it once on the full matrix
%   introduces no leakage -- each subject's value depends only on that
%   subject's own row. Covariates are still regressed out of that feature per
%   fold, train-only, so the baseline is controlled the same way the model is.
%
%   See also foldPreprocess, residualizeFold, PLSR_neuroimaging_pipeline.

n = size(X,1);

if isa(opts.globalFun,'function_handle')
    gf = opts.globalFun(X);
else
    switch lower(opts.globalFun)
        case 'median'
            gf = median(X,2);
        otherwise
            gf = mean(X,2);
    end
end
gf = gf(:);

isClass = strcmpi(kind,'classify');

vals   = nan(opts.nRepeats, opts.outerK);
valsPR = nan(opts.nRepeats, opts.outerK);

for r = 1:opts.nRepeats

    if isClass
        try
            cvO = cvpartition(Y,'KFold',opts.outerK,'Stratify',true);
        catch
            cvO = cvpartition(n,'KFold',opts.outerK);
        end
    else
        cvO = cvpartition(n,'KFold',opts.outerK);
    end

    for k = 1:opts.outerK

        tr = training(cvO,k);
        te = test(cvO,k);

        ytr = Y(tr);
        yte = Y(te);

        if isClass && (numel(unique(ytr)) < 2 || numel(unique(yte)) < 2)
            continue
        end
        if ~isClass && (numel(ytr) < 3 || numel(yte) < 2)
            continue
        end

        Ctr = []; Cte = [];
        if ~isempty(C), Ctr = C(tr,:); Cte = C(te,:); end

        [gtr, gte] = residualizeFold(gf(tr), gf(te), Ctr, Cte);

        if std(gtr) == 0
            continue
        end

        try
            if isClass
                mdl = fitglm(gtr, ytr, 'Distribution','binomial');
                sc  = predict(mdl, gte);
                if numel(unique(sc)) < 2 || ~all(isfinite(sc))
                    continue
                end
                [~,~,~,vals(r,k)]   = perfcurve(yte, sc, 1);
                [~,~,~,valsPR(r,k)] = perfcurve(yte, sc, 1, 'xCrit','reca','yCrit','prec');
            else
                mdl  = fitlm(gtr, ytr);
                yhat = predict(mdl, gte);
                d    = sum((yte - mean(ytr)).^2);
                if d <= 0
                    continue
                end
                vals(r,k) = 1 - sum((yte - yhat).^2) / d;
            end
        catch
            % leave this fold as NaN
        end

    end
end

if isClass
    out.AUC    = mean(vals(:),  'omitnan');
    out.AUC_PR = mean(valsPR(:),'omitnan');
else
    out.Q2 = mean(vals(:),'omitnan');
end

end
