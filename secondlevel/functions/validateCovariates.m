function [C, covInfo] = validateCovariates(opts, n, pipelineName)
% validateCovariates  Resolve and check opts.covariates once, up front.
%
%   [C, covInfo] = validateCovariates(opts, n, pipelineName)
%
%   Reads the optional opts.covariates field, validates it against the sample
%   size, and returns both the matrix and a provenance struct for the results.
%   Called once at the top of each pipeline; the resulting C is then passed
%   EXPLICITLY to every helper alongside X, never re-read from opts.
%
%   Inputs
%     opts          options struct; reads .covariates and .covariateNames
%     n             number of subjects (rows of X)
%     pipelineName  char, used in error and warning messages
%
%   Outputs
%     C        [n x c] numeric covariate matrix, or [] when none supplied
%     covInfo  struct with .used, .names, .nCov, .rank
%
%   Validation performed
%     - numeric and 2-D
%     - exactly n rows (a transposed matrix is the common mistake, so that case
%       gets its own message)
%     - all values finite; missing data must be handled upstream, consistent
%       with the rest of the family's input contract
%     - warns on columns that are constant across the whole sample (they carry
%       no information and are absorbed by the intercept)
%     - warns when the covariate design is rank deficient
%     - warns when the covariate count is large relative to n, since each
%       covariate costs a degree of freedom in every training fold
%
%   Do NOT supply a column of ones. The intercept is handled internally by
%   residualizeFold, which centers covariates on the training-fold mean.
%
%   See also residualizeFold, foldPreprocess.

C = [];
covInfo = struct('used',false,'names',{{}},'nCov',0,'rank',0);

if ~isfield(opts,'covariates') || isempty(opts.covariates)
    return
end

C = opts.covariates;

if ~isnumeric(C) || ~ismatrix(C)
    error([pipelineName ':covariatesType'], ...
        'opts.covariates must be a numeric [n x nCov] matrix. Encode categorical covariates as dummy columns before calling.');
end

if size(C,1) ~= n
    if size(C,2) == n
        error([pipelineName ':covariatesTransposed'], ...
            'opts.covariates is [%d x %d] but X has %d rows. It looks transposed; expected one ROW per subject.', ...
            size(C,1), size(C,2), n);
    end
    error([pipelineName ':covariatesRows'], ...
        'opts.covariates has %d rows but X has %d. One row per subject is required.', size(C,1), n);
end

if ~all(isfinite(C(:)))
    error([pipelineName ':covariatesNonFinite'], ...
        'opts.covariates contains NaN or Inf. Handle missing values upstream (impute or drop subjects).');
end

nCov = size(C,2);

if isfield(opts,'covariateNames') && ~isempty(opts.covariateNames)
    names = cellstr(opts.covariateNames);
    if numel(names) ~= nCov
        error([pipelineName ':covariateNames'], ...
            'opts.covariateNames has %d entries but opts.covariates has %d columns.', numel(names), nCov);
    end
else
    names = arrayfun(@(i) sprintf('cov%d',i), 1:nCov, 'UniformOutput', false);
end

constCols = std(C,0,1) == 0;
if any(constCols)
    warning([pipelineName ':covariateConstant'], ...
        'Covariate(s) %s are constant across all subjects and carry no information; they will be dropped per fold.', ...
        strjoin(names(constCols), ', '));
end

r = rank([ones(n,1) C]);
if r < nCov + 1
    warning([pipelineName ':covariateRankDeficient'], ...
        'The covariate design [1 C] has rank %d but %d columns: some covariates are collinear. The nuisance fit uses a minimum-norm solution, so individual coefficients are not uniquely identified (the residuals still are).', ...
        r, nCov + 1);
end

if nCov > n/5
    warning([pipelineName ':covariateCount'], ...
        'Using %d covariates with only %d subjects. Each covariate costs a degree of freedom in every training fold; consider reducing the nuisance model.', ...
        nCov, n);
end

covInfo = struct('used',true,'names',{names},'nCov',nCov,'rank',r);

end
