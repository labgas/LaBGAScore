function lambdaGrid = enetLambdaGrid(opts)
% enetLambdaGrid  Resolve the Elastic Net lambda path from an options struct.
%
%   lambdaGrid = enetLambdaGrid(opts) returns opts.lambdaGrid if present and
%   non-empty, otherwise the family default logspace(-3,1,25).
%
%   This exists so that the point estimate (ENet_neuroimaging_pipeline), the
%   bootstrap CI (bootstrapOOB_ENet) and the permutation null (quickCV_ENet,
%   quickCV_ENet_PR) all search the SAME lambda path. They previously did not:
%   the main nested CV honoured opts.lambdaGrid while the CI and null hardcoded
%   'NumLambda',25 / 'LambdaRatio',1e-3, so the three quantities reported side
%   by side in the results struct were computed over different grids.
%
%   Input
%     opts  struct; only the optional field .lambdaGrid is read.
%
%   Output
%     lambdaGrid  [1 x nLambda] positive lambda values.
%
%   Note on ordering: lassoglm preserves the order of a user-supplied Lambda
%   vector, so column l of B corresponds to lambdaGrid(l). Do not assume the
%   descending-order convention lassoglm uses when it generates its own path.
%
%   See also lassoglm, ENet_neuroimaging_pipeline, quickCV_ENet, bootstrapOOB_ENet.

if isfield(opts,'lambdaGrid') && ~isempty(opts.lambdaGrid)
    lambdaGrid = opts.lambdaGrid;
else
    lambdaGrid = logspace(-3,1,25);
end

end
