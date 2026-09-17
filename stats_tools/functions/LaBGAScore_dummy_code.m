function [X, levels, catvec] = LaBGAScore_dummy_code(v, varargin)
% Dummy-code a phenotype column into k-1 indicator columns.
%
% *USAGE*
%
% X                    = LaBGAScore_dummy_code(v)
% [X, levels, catvec]  = LaBGAScore_dummy_code(v, 'reference', 'KUL')
%
% *WHY*
%
% A k-level factor needs k-1 columns. Coding three sites as a single -1/0/1
% column - as proj_discoverie's model_3a does for 'center' - treats them as
% ordered and spends only one degree of freedom, so it removes only part of
% the between-site variance and imposes an ordering that does not exist. Any
% nuisance adjustment for an unordered factor should go through this.
%
% Accepts numeric, logical, char, cellstr, string or categorical input and
% returns double indicator columns, dropping one level as the reference.
%
% *INPUTS*
%
%   v            column of values, one per subject
%
% *OPTIONAL INPUTS*
%
%   'reference'  level to drop. Default: the first level, which for a
%                categorical is its first category and otherwise sorted order.
%
% *OUTPUTS*
%
%   X            n x (k-1) double, indicator columns for levels 2..k
%   levels       the k level names, in order; levels{1} is the reference
%   catvec       the input as a categorical, useful for crosstab
%
% *SEE ALSO*
%
% LaBGAScore_decoding_SVM_between_subjects, prep_3a_run_second_level_regression_and_save
%
% -------------------------------------------------------------------------
% Lukas Van Oudenhove, KU Leuven, September 2026
% -------------------------------------------------------------------------

p = inputParser;
p.addParameter('reference', '', @(x) ischar(x) || isstring(x) || isnumeric(x));
p.parse(varargin{:});

if iscategorical(v)
    catvec = removecats(v(:));
elseif isnumeric(v) || islogical(v)
    catvec = categorical(v(:));
else
    catvec = categorical(cellstr(string(v(:))));
end

levels = categories(catvec);
k = numel(levels);

if k < 2
    X = zeros(numel(catvec), 0);
    return
end

ref = p.Results.reference;
if ~isempty(ref)
    ref = char(string(ref));
    if ~ismember(ref, levels)
        error('LaBGAScore_dummy_code:badref', ...
            'reference ''%s'' is not a level. Levels: %s', ref, strjoin(levels', ', '));
    end
    levels = [{ref}; levels(~strcmp(levels, ref))];
end

X = zeros(numel(catvec), k-1);
for i = 2:k
    X(:, i-1) = double(catvec == levels{i});
end

end
