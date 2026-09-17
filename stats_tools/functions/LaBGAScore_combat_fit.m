function est = LaBGAScore_combat_fit(dat, batch, varargin)
% Fit ComBat harmonization parameters on a TRAINING set only
%
% :Usage:
% ::
%     est = LaBGAScore_combat_fit(dat, batch)
%     est = LaBGAScore_combat_fit(dat, batch, 'parametric', 1, 'ref', refLabel)
%
% :Inputs:
%   **dat:**      p x n data matrix (features x subjects), training set only
%   **batch:**    n x 1 batch/site labels (numeric, cellstr or categorical)
%
% :Optional inputs:
%   **'parametric':** 1 (default) or 0, empirical Bayes adjustment type
%   **'ref':**        label of the reference batch; [] harmonizes to the grand mean
%
% :Outputs:
%   **est:** struct of fitted parameters for LaBGAScore_combat_apply:
%            .batch_labels, .grand_mean, .var_pooled, .gamma_star, .delta_star, .ref
%
% :Notes:
% Written for CROSS-VALIDATED DECODING, where ComBat must be fitted on the
% training fold and applied to held-out subjects. It therefore fits with NO
% covariates (mod = []), deliberately, and rejects a non-empty mod.
%
% The reason is structural, not a simplification. combat.m builds
%     stand_mean = grand_mean + (design with batch columns zeroed) * B_hat
% so with covariates the standardization of a subject depends on that
% subject's covariate values. Applying fitted parameters to a held-out subject
% would then require that subject's covariates - and in decoding the covariate
% of interest is the outcome being predicted, which is exactly what must not be
% used. The R reference implementation takes the same position:
% neuroCombatFromTraining(dat, batch, estimates, mod = NULL).
%
% Consequence worth stating: with mod = [], ComBat removes each batch's mean
% and variance. A batch whose subjects are all one class therefore has its
% class signal removed along with its site effect, because the two are not
% separable within that batch. Do not use this on such a batch and expect the
% class to remain decodable there.
%
% The empirical Bayes estimation itself is delegated to combat.m (Jfortin1/
% ComBatHarmonization) rather than reimplemented. Only grand_mean and
% var_pooled, which combat.m does not return, are recomputed here, using the
% same formulas (combat.m lines 136-152). LaBGAScore_combat_apply on the
% training set reproduces combat.m's own output exactly; that identity is the
% correctness test for this pair of functions.
%
% :See also: LaBGAScore_combat_apply, combat
%
% -------------------------------------------------------------------------
% Author: Lukas Van Oudenhove
% Date: September, 2026
% -------------------------------------------------------------------------
% LaBGAScore_combat_fit.m         v1.0
% -------------------------------------------------------------------------

% -------------------------------- options --------------------------------

parametric = 1;
ref        = [];
modarg     = [];

for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'parametric', parametric = double(logical(varargin{i+1}));
        case 'ref',        ref = varargin{i+1};
        case 'mod',        modarg = varargin{i+1};
        otherwise, error('Unrecognized option ''%s''.', varargin{i});
    end
end

if ~isempty(modarg)
    error(['%s'], sprintf(['\nLaBGAScore_combat_fit does not accept a non-empty mod.\n' ...
        'Fitted parameters could then only be applied to a subject whose covariates are\n' ...
        'known, and in decoding the covariate of interest is the outcome. See the notes\n' ...
        'in the help text.\n']));
end

if isempty(which('combat'))
    error('combat.m is not on the path. Add ComBatHarmonization/Matlab/scripts.');
end

% ---------------------------- batch bookkeeping --------------------------

batch = batch(:);
labels_all = cellstr(string(batch));
[batch_labels, ~, batch_idx] = unique(labels_all, 'stable');
batch_idx = double(batch_idx);
n_batch   = numel(batch_labels);
n_array   = numel(batch_idx);

if size(dat,2) ~= n_array
    error('dat has %d columns but batch has %d entries.', size(dat,2), n_array);
end

ref_code = [];
if ~isempty(ref)
    ref_code = find(strcmp(batch_labels, char(string(ref))));
    if isempty(ref_code)
        error('reference batch ''%s'' is not present among the training batches (%s).', ...
            char(string(ref)), strjoin(batch_labels(:)', ', '));
    end
end

% ----------------- empirical Bayes parameters, from combat.m -------------

cb_args = {double(dat), batch_idx, [], parametric};
if ~isempty(ref_code), cb_args = [cb_args {'ref', ref_code}]; end
[~, gamma_star, delta_star] = combat(cb_args{:});

% -------- standardization constants, recomputed as in combat.m -----------
% combat.m does not return these, so they are rebuilt here with the same
% formulas. design is batch dummies only, because mod is empty by contract.

batchmod = dummyvar({categorical(batch_idx)});
design   = batchmod;
if ~isempty(ref_code)
    design(:, ref_code) = 1;      % combat.m line 114
end

n_batches = accumarray(batch_idx, 1)';
B_hat     = (design' * design) \ (design' * double(dat)');

if isempty(ref_code)
    grand_mean = (n_batches / n_array) * B_hat(1:n_batch, :);
    var_pooled = ((double(dat) - (design * B_hat)').^2) * repmat(1/n_array, n_array, 1);
else
    grand_mean = B_hat(ref_code, :);
    ref_sel    = batch_idx == ref_code;
    ref_dat    = double(dat(:, ref_sel));
    ref_n      = n_batches(ref_code);
    var_pooled = ((ref_dat - (design(ref_sel,:) * B_hat)').^2) * repmat(1/ref_n, ref_n, 1);
end

% --------------------------------- output --------------------------------

est = struct();
est.batch_labels = batch_labels;
est.grand_mean   = grand_mean;     % 1 x p
est.var_pooled   = var_pooled;     % p x 1
est.gamma_star   = gamma_star;     % n_batch x p
est.delta_star   = delta_star;     % n_batch x p
est.ref          = ref;
est.ref_code     = ref_code;
est.parametric   = parametric;
est.n_per_batch  = n_batches;

end
