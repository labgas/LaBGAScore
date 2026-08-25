function [beta, intercept, alphaSel, lambdaSel, tuneInfo] = ...
    selectENetHyperparams(XtrRaw, ytr, Ctr, cvIn, alphaGrid, lambdaGrid, rule, scaleMode)
% selectENetHyperparams  Fold-averaged Elastic Net hyperparameter selection.
%
%   [beta, intercept, alphaSel, lambdaSel, tuneInfo] = selectENetHyperparams( ...
%       XtrRaw, ytr, Ctr, cvIn, alphaGrid, lambdaGrid, rule, scaleMode)
%
%   Selects (alpha, lambda) by inner cross-validation on the TRAINING FOLD ONLY,
%   averaging AUC across inner folds before comparing candidates, then refits
%   once on the full training fold at the winning pair.
%
%   Inputs
%     XtrRaw      [nTr x p]  RAW training-fold features (not yet preprocessed)
%     ytr         [nTr x 1]  binary 0/1 labels
%     Ctr         [nTr x c]  training-fold covariates, or []
%     cvIn        cvpartition over the training fold
%     alphaGrid   [1 x nA]   Elastic Net mixing values
%     lambdaGrid  [1 x nL]   regularization path (see enetLambdaGrid)
%     rule        '1se' (default) | 'max'
%     scaleMode   'zscore' | 'center' | 'none'
%
%   Outputs
%     beta        [p x 1]  coefficients from the single refit, on the scale of
%                          foldPreprocess(XtrRaw, ., Ctr, ., scaleMode)
%     intercept   scalar
%     alphaSel    selected alpha (NaN if no inner fold produced a usable model)
%     lambdaSel   selected lambda
%     tuneInfo    struct with .meanAUC [nA x nL], .seAUC, .nValidFolds, .rule
%
%   Raw features go in, not preprocessed ones, because each inner fold performs
%   its OWN train-only preprocessing. Reusing the outer fold's scaling and
%   nuisance coefficients would leak outer-training statistics into every inner
%   validation fold -- harmless-looking for z-scoring, but exactly the leak that
%   fold-wise covariate regression exists to remove, one level down.
%
%   Why the selection rule changed
%     The previous inline search in ENet_neuroimaging_pipeline compared a running
%     maximum against the AUC of a SINGLE inner validation fold, over a flattened
%     alpha x fold x lambda grid, with the running best never reset per alpha.
%     Hyperparameters were therefore chosen by the luckiest fold. Measured on
%     true-null data with the production grid, that rule picked a median lambda
%     of 0.068 against 1.47 for fold-averaged selection -- roughly twenty times
%     too little regularization -- and gave a noticeably more variable outer AUC
%     (sd 0.208 vs 0.160). It did NOT inflate the outer AUC itself, because the
%     outer test fold is held out from selection; the damage is unstable,
%     under-regularized models and hence unstable weight maps.
%     It also refit the whole training fold inside the innermost loop on every
%     new running maximum, instead of once at the end.
%
%   See also lassoglm, enetLambdaGrid, foldPreprocess, ENet_neuroimaging_pipeline.

if nargin < 7 || isempty(rule),      rule = '1se';         end
if nargin < 8 || isempty(scaleMode), scaleMode = 'zscore'; end

nA = numel(alphaGrid);
nL = numel(lambdaGrid);
K  = cvIn.NumTestSets;

foldAUC = nan(nA, K, nL);

for f = 1:K

    tr = training(cvIn,f);
    va = test(cvIn,f);

    y2 = ytr(tr);
    yv = ytr(va);

    if numel(unique(y2)) < 2 || numel(unique(yv)) < 2
        continue
    end

    C2 = []; Cv = [];
    if ~isempty(Ctr), C2 = Ctr(tr,:); Cv = Ctr(va,:); end

    % train-only preprocessing for THIS inner fold
    [X2, Xv] = foldPreprocess(XtrRaw(tr,:), XtrRaw(va,:), C2, Cv, scaleMode);

    for a = 1:nA
        try
            [B,FI] = lassoglm( ...
                X2, y2, 'binomial', ...
                'Alpha', alphaGrid(a), ...
                'Lambda', lambdaGrid, ...
                'Standardize', false);
        catch
            continue
        end

        % lassoglm preserves the supplied Lambda order, so column l <-> lambdaGrid(l)
        for l = 1:size(B,2)
            scores = Xv*B(:,l) + FI.Intercept(l);
            if ~all(isfinite(scores)) || numel(unique(scores)) < 2
                continue
            end
            try
                [~,~,~,foldAUC(a,f,l)] = perfcurve(yv, scores, 1);
            catch
                % leave as NaN
            end
        end
    end

end

% ---- average across inner folds BEFORE comparing candidates ----
meanAUC = reshape(mean(foldAUC, 2, 'omitnan'),   [nA nL]);
nValid  = reshape(sum(isfinite(foldAUC), 2),     [nA nL]);
sdAUC   = reshape(std(foldAUC, [], 2, 'omitnan'),[nA nL]);
seAUC   = sdAUC ./ sqrt(max(nValid,1));

tuneInfo = struct('meanAUC',meanAUC,'seAUC',seAUC,'nValidFolds',nValid,'rule',rule);

p = size(XtrRaw,2);

if all(~isfinite(meanAUC(:)))
    beta      = nan(p,1);
    intercept = NaN;
    alphaSel  = NaN;
    lambdaSel = NaN;
    return
end

[bestVal, bestIdx] = max(meanAUC(:));

switch lower(rule)
    case 'max'
        selIdx = bestIdx;

    otherwise % '1se'
        thresh = bestVal - seAUC(bestIdx);
        if ~isfinite(thresh)
            thresh = bestVal;
        end
        cand = find(meanAUC >= thresh);
        [~, lCand] = ind2sub([nA nL], cand);
        % most regularized = largest lambda; break ties on the best mean AUC
        maxLambdaIdx = max(lCand);
        tied = cand(lCand == maxLambdaIdx);
        [~,w] = max(meanAUC(tied));
        selIdx = tied(w);
end

[ia, il] = ind2sub([nA nL], selIdx);
alphaSel  = alphaGrid(ia);
lambdaSel = lambdaGrid(il);

% ---- single refit on the full training fold ----
% Preprocessing the training fold against itself reproduces exactly the XtrP the
% caller obtained from foldPreprocess(XtrRaw, Xte, Ctr, Cte, scaleMode), since
% every constant is estimated from the training rows alone.
XtrP = foldPreprocess(XtrRaw, [], Ctr, [], scaleMode);

try
    [Bf, FIf] = lassoglm( ...
        XtrP, ytr, 'binomial', ...
        'Alpha', alphaSel, ...
        'Lambda', lambdaSel, ...
        'Standardize', false);
    beta      = Bf(:,1);
    intercept = FIf.Intercept(1);
catch
    beta      = nan(p,1);
    intercept = NaN;
end

end
