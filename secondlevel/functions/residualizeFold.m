function [XtrR, XteR, nuis] = residualizeFold(Xtr, Xte, Ctr, Cte)
% residualizeFold  Train-only nuisance regression for one cross-validation fold.
%
%   [XtrR, XteR, nuis] = residualizeFold(Xtr, Xte, Ctr, Cte)
%
%   Regresses the nuisance covariates out of the feature matrix, estimating the
%   nuisance coefficients on the TRAINING fold only and applying them to both
%   the training and the test fold. This mirrors the contract of applyScaling:
%   every parameter is fit on train and merely applied to test, so no test-fold
%   information reaches the model.
%
%   This is the leakage-free replacement for residualizing the full feature
%   matrix before cross-validation. Fitting the nuisance model on all subjects
%   lets each test fold influence the nuisance coefficients that are later used
%   to transform it, which inflates every downstream performance estimate.
%
%   Inputs
%     Xtr   [nTr x p]  training features (raw, not yet scaled)
%     Xte   [nTe x p]  test features
%     Ctr   [nTr x c]  training covariates, or [] to disable
%     Cte   [nTe x c]  test covariates, or []
%
%   Outputs
%     XtrR  [nTr x p]  training residuals
%     XteR  [nTe x p]  test residuals, using TRAIN means and TRAIN coefficients
%     nuis  struct describing the fit:
%             .used        false when Ctr was empty (pure passthrough)
%             .muC         training-fold covariate means
%             .keep        columns retained (non-constant within this fold)
%             .Bc          [nKeep x p] nuisance coefficients
%             .rank        effective rank of the centered covariate matrix
%             .varRemoved  [1 x p] fraction of per-feature variance removed
%
%   Implementation notes
%     - Covariates are centered on the TRAINING mean, and only the covariate-
%       attributable part is subtracted. The intercept is therefore absorbed
%       rather than removed, so feature means survive for applyScaling to
%       handle and residualization is an exact no-op when the covariates carry
%       no variance. This is what makes the empty-covariate path bit-identical
%       to the previous behaviour.
%     - A truncated SVD (not backslash, not blind pinv) supplies the minimum-
%       norm solution under collinearity, without a warning and without an
%       arbitrary choice of which collinear column to keep.
%     - Columns that are constant WITHIN this fold are dropped. A covariate can
%       be non-constant overall yet constant inside a small CV fold.
%
%   See also foldPreprocess, applyScaling, validateCovariates.

if isempty(Ctr)
    XtrR = Xtr;
    XteR = Xte;
    nuis = struct('used',false,'muC',[],'keep',[],'Bc',[],'rank',0,'varRemoved',[]);
    return
end

if size(Ctr,1) ~= size(Xtr,1)
    error('residualizeFold:rowMismatchTrain', ...
        'Ctr has %d rows but Xtr has %d. Covariates must be sliced with the same index as X.', ...
        size(Ctr,1), size(Xtr,1));
end
if ~isempty(Xte) && size(Cte,1) ~= size(Xte,1)
    error('residualizeFold:rowMismatchTest', ...
        'Cte has %d rows but Xte has %d. Covariates must be sliced with the same index as X.', ...
        size(Cte,1), size(Xte,1));
end

muC = mean(Ctr,1);
Cc  = Ctr - muC;

keep = std(Cc,0,1) > 0;
Cc   = Cc(:,keep);

if isempty(Cc)
    XtrR = Xtr;
    XteR = Xte;
    nuis = struct('used',false,'muC',muC,'keep',keep,'Bc',[],'rank',0,'varRemoved',zeros(1,size(Xtr,2)));
    return
end

[U,S,V] = svd(Cc,'econ');
sv  = diag(S);
tol = max(size(Cc)) * eps(sv(1));
k   = sum(sv > tol);

U  = U(:,1:k);
Sk = diag(sv(1:k));
Vk = V(:,1:k);

P    = U' * Xtr;              % [k x p] projection onto the covariate subspace
XtrR = Xtr - U*P;             % better conditioned than Cc*Bc, numerically identical
Bc   = Vk * (Sk \ P);         % [nKeep x p] coefficients in covariate space

if isempty(Xte)
    XteR = Xte;
else
    XteR = Xte - (Cte(:,keep) - muC(keep)) * Bc;
end

nuis = struct( ...
    'used',       true, ...
    'muC',        muC, ...
    'keep',       keep, ...
    'Bc',         Bc, ...
    'rank',       k, ...
    'varRemoved', 1 - var(XtrR,0,1) ./ max(var(Xtr,0,1), eps));

end
