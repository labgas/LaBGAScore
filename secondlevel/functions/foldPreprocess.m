function [XtrP, XteP, pre] = foldPreprocess(Xtr, Xte, Ctr, Cte, scaleMode)
% foldPreprocess  Per-fold, train-only preprocessing for the ML pipeline family.
%
%   [XtrP, XteP, pre] = foldPreprocess(Xtr, Xte, Ctr, Cte, scaleMode)
%
%   Composes the two train-only transformations every pipeline fold needs, in
%   the one order the family agrees on:
%
%       1. nuisance regression  (residualizeFold)  -- covariates out
%       2. scaling              (applyScaling)     -- z-score / center / none
%
%   Both steps estimate their parameters on the TRAINING fold and merely apply
%   them to the test fold.
%
%   Inputs
%     Xtr, Xte   feature matrices for this fold (raw)
%     Ctr, Cte   covariates for this fold, or [] for none
%     scaleMode  'zscore' | 'center' | 'none'
%
%   Outputs
%     XtrP, XteP  preprocessed features
%     pre         struct with .nuis (from residualizeFold) and .sc (scaling)
%
%   Why the order is residualize-then-scale
%     For the TRAINING fold the two orders are algebraically equivalent, because
%     the nuisance projection leaves column means intact. They differ in three
%     ways that matter:
%       - The test fold genuinely differs, since the scaling constants differ.
%       - The scale constants are computed on the residualized data, so the
%         model's inputs really do have unit variance in the space it operates
%         in. This matters for ENet, which calls lassoglm with
%         'Standardize',false: the L1/L2 penalty applies directly to the
%         supplied columns, so mis-scaled inputs silently reweight the penalty
%         per feature.
%       - Residualization lowers rank(Xtr) by up to rank(C), so capLV must see
%         the post-residualization matrix or it will request more latent
%         variables than the data support. Residualizing first gets that free.
%     It also matches the standard neuroimaging order: confound regression, then
%     normalization.
%
%   Why this is one function rather than two calls at each site
%     The order-of-operations decision lives here and nowhere else, and no call
%     site can silently omit the residualization step. With Ctr = Cte = [] this
%     degenerates exactly to applyScaling, which is what keeps every existing
%     3-argument pipeline call bit-identical to its previous behaviour.
%
%   IMPORTANT: covariates are always passed explicitly alongside X, never read
%   from an opts struct inside a helper. Several pipeline stages subset subjects
%   (bootstrap resampling, learning-curve subsampling), and a full-length
%   covariate matrix reaching a subsetted X would be silently misaligned rather
%   than an error. The row check below turns that class of bug into a hard fail.
%
%   See also residualizeFold, applyScaling, validateCovariates.

if nargin < 5 || isempty(scaleMode)
    scaleMode = 'zscore';
end

[XtrP, XteP, pre.nuis] = residualizeFold(Xtr, Xte, Ctr, Cte);
[XtrP, XteP, pre.sc]   = applyScaling(XtrP, XteP, scaleMode);

end
