function out = LaBGAScore_combat_apply(dat, batch, est)
% Apply ComBat parameters fitted on a training set to held-out data
%
% :Usage:
% ::
%     out = LaBGAScore_combat_apply(dat, batch, est)
%
% :Inputs:
%   **dat:**    p x n data matrix (features x subjects) to harmonize
%   **batch:**  n x 1 batch/site labels for those subjects
%   **est:**    struct returned by LaBGAScore_combat_fit
%
% :Outputs:
%   **out:**    p x n harmonized data
%
% :Notes:
% No outcome or covariate information is used, by construction - only the
% subject's batch. That is what makes this safe to run on held-out subjects
% inside a cross-validation loop.
%
% Every batch present here must have been present when est was fitted: a batch
% the training set never saw has no estimated parameters, so leave-one-site-out
% designs are NOT supported. The R reference implementation raises the same
% error (neuroCombat.R, 'missing.levels').
%
% The transform inverts combat.m's own steps (lines 160, 226, 229) with the
% training-set constants:
%     s   = (dat - grand_mean) ./ sqrt(var_pooled)
%     s   = (s - gamma_star(b,:)) ./ sqrt(delta_star(b,:))
%     out = s .* sqrt(var_pooled) + grand_mean
% Applying this to the training set itself reproduces combat.m's output
% exactly; that identity is the correctness test for this pair of functions.
%
% :See also: LaBGAScore_combat_fit, combat
%
% -------------------------------------------------------------------------
% Author: Lukas Van Oudenhove
% Date: September, 2026
% -------------------------------------------------------------------------
% LaBGAScore_combat_apply.m         v1.0
% -------------------------------------------------------------------------

batch = batch(:);
labels_here = cellstr(string(batch));

if size(dat,2) ~= numel(labels_here)
    error('dat has %d columns but batch has %d entries.', size(dat,2), numel(labels_here));
end
if size(dat,1) ~= numel(est.grand_mean)
    error('dat has %d features but est was fitted on %d.', size(dat,1), numel(est.grand_mean));
end

missing = setdiff(unique(labels_here), est.batch_labels);
if ~isempty(missing)
    error(['%s'], sprintf(['\nBatch(es) {%s} are not in the fitted estimates {%s}.\n' ...
        'ComBat cannot harmonize a site it never saw during fitting, so leave-one-site-out\n' ...
        'cross-validation is not supported. Stratify the folds by site instead.\n'], ...
        strjoin(missing(:)', ', '), strjoin(est.batch_labels(:)', ', ')));
end

sd_pooled  = sqrt(est.var_pooled);              % p x 1
stand_mean = est.grand_mean(:);                 % p x 1
out        = double(dat);

for k = 1:numel(est.batch_labels)
    sel = strcmp(labels_here, est.batch_labels{k});
    if ~any(sel), continue; end
    s = (out(:,sel) - stand_mean) ./ sd_pooled;
    s = (s - est.gamma_star(k,:)') ./ sqrt(est.delta_star(k,:)');
    out(:,sel) = s .* sd_pooled + stand_mean;
end

end
