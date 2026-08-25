function [ytrR, yteR] = residualizeY(ytr, yte, Ctr, Cte, doIt)
% residualizeY  Optional train-only nuisance regression of a continuous outcome.
%
%   [ytrR, yteR] = residualizeY(ytr, yte, Ctr, Cte, doIt)
%
%   Companion to residualizeFold for the outcome side. Only meaningful for a
%   CONTINUOUS outcome, so only PLSR_neuroimaging_pipeline uses it: the three
%   classification pipelines binarize Y as double(Y == max(Y)), and a
%   residualized Y is no longer binary, which would break perfcurve, lassoglm
%   with 'binomial', and the scores > 0.5 decision rule.
%
%   Inputs
%     ytr, yte  training and test outcome for this fold
%     Ctr, Cte  fold covariates, or []
%     doIt      logical; false (the default in PLSR) returns the inputs unchanged
%
%   Outputs
%     ytrR, yteR  residualized outcome, using TRAIN means and TRAIN coefficients
%
%   Interpretation
%     With doIt = false and covariates supplied, X is confound-free but Y is
%     not, so the nuisance-explained part of Y is structurally unexplainable and
%     Q2 is deflated relative to a partial-Q2. With doIt = true both sides are
%     residualized and Q2 becomes a partial Q2: variance in Y explained by X
%     after removing the covariates from both, the regression analogue of a
%     partial correlation. The two are not comparable, so record which was used.
%
%     The Q2 denominator downstream uses mean(ytrain) of whatever this returns,
%     which stays correct in both cases because the training mean is recomputed
%     from the returned values.
%
%   See also residualizeFold, foldPreprocess, PLSR_neuroimaging_pipeline.

if nargin < 5 || isempty(doIt) || ~doIt || isempty(Ctr)
    ytrR = ytr;
    yteR = yte;
    return
end

[ytrR, yteR] = residualizeFold(ytr(:), yte(:), Ctr, Cte);

end
