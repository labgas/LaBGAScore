function AUC = quickCV_PLSDA(X, Y, C, opts)
% quickCV_PLSDA  Fast tuned cross-validated ROC AUC for PLS-DA.
%
%   AUC = quickCV_PLSDA(X, Y, C, opts)
%
%   Thin wrapper over quickCV_PLSDA_core with metric 'roc'. Used by the
%   permutation, learning-curve and matched-observed stages of
%   PLSDA_neuroimaging_pipeline.
%
%   Pass C = [] when no covariates are in use. C must be sliced with the same
%   subject index as X; see foldPreprocess for why covariates travel as an
%   explicit argument rather than inside opts.
%
%   See also quickCV_PLSDA_core, quickCV_PLSDA_PR, PLSDA_neuroimaging_pipeline.

if nargin < 4
    error('quickCV_PLSDA:nargin', ...
        'quickCV_PLSDA now requires (X, Y, C, opts). Pass C = [] when no covariates are used.');
end

AUC = quickCV_PLSDA_core(X, Y, C, opts, 'roc');

end
