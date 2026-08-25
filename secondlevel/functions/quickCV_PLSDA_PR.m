function AUC_PR = quickCV_PLSDA_PR(X, Y, C, opts)
% quickCV_PLSDA_PR  Fast tuned cross-validated precision-recall AUC for PLS-DA.
%
%   AUC_PR = quickCV_PLSDA_PR(X, Y, C, opts)
%
%   Thin wrapper over quickCV_PLSDA_core with metric 'pr'. Note that the chance
%   level for precision-recall AUC is the positive-class prevalence, not 0.5.
%
%   Pass C = [] when no covariates are in use.
%
%   See also quickCV_PLSDA_core, quickCV_PLSDA, PLSDA_neuroimaging_pipeline.

if nargin < 4
    error('quickCV_PLSDA_PR:nargin', ...
        'quickCV_PLSDA_PR now requires (X, Y, C, opts). Pass C = [] when no covariates are used.');
end

AUC_PR = quickCV_PLSDA_core(X, Y, C, opts, 'pr');

end
