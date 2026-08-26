function [dat, lo, hi] = maskToSignificant(stat_vec, sig)
% maskToSignificant  Prepare a statistic vector so a range threshold keeps only
% the significant voxels.
%
%   [dat, lo, hi] = maskToSignificant(stat_vec, sig)
%
%   Returns stat_vec with the NON-significant voxels moved to a sentinel value
%   below lo, together with a window [lo hi] that brackets exactly the
%   significant values. Thresholding the result with
%
%       obj.threshold([lo hi], 'raw-between', 'k', k)
%
%   then keeps precisely the significant voxels while still applying the
%   cluster-extent rule.
%
%   Inputs
%     stat_vec  [nVox x 1] statistic values for every voxel
%     sig       [nVox x 1] logical mask of significant voxels
%
%   Outputs
%     dat       stat_vec with non-significant voxels set to the sentinel
%     lo, hi    threshold window bracketing the significant values
%
%   Why this is needed
%     The previous code thresholded on
%         [min(stat_vec(sig)) - 0.001, max(stat_vec(sig)) + 0.001]
%     applied to the UNMASKED statistic vector. 'raw-between' selects on value,
%     so that retains every voxel whose statistic happens to fall inside the
%     window, significant or not. It is equivalent to masking only when the
%     p-value is monotonic in the statistic. For voxelwise permutation
%     p-values it is not: each voxel is compared against its own null, so a
%     non-significant voxel can easily carry a statistic between the smallest
%     and largest significant ones, and it was being kept.
%
%     The 0.001 pad had the same flavour of problem: a hardcoded absolute
%     tolerance applied to statistics whose scale is arbitrary (AUC, TFCE,
%     t-values). It is replaced here by a relative tolerance.
%
%   The sentinel is used rather than zero because significant values may
%   straddle zero for a signed statistic, in which case zeroed voxels would sit
%   inside the retained window.
%
%   See also thresholded_fmri_data_from_statistic_image,
%   thresholded_fmri_data_from_pval_nii.

stat_vec = stat_vec(:);
sig      = logical(sig(:));

if numel(sig) ~= numel(stat_vec)
    error('maskToSignificant:sizeMismatch', ...
        'stat_vec has %d elements but sig has %d.', numel(stat_vec), numel(sig));
end

vals = stat_vec(sig);
vals = vals(isfinite(vals));

if isempty(vals)
    error('maskToSignificant:noFiniteSignificantValues', ...
        'No finite statistic values among the significant voxels.');
end

lo0  = min(vals);
hi0  = max(vals);
span = hi0 - lo0;

scale = max(abs([lo0 hi0]));
if ~isfinite(scale) || scale == 0
    scale = 1;
end

% relative tolerance, so the window does not depend on the statistic's units
tol = max(span, scale) * 1e-6;

lo = lo0 - tol;
hi = hi0 + tol;

% strictly below lo, so the sentinel can never fall inside the window
sentinel = lo - max(span, scale) - 1;

dat = stat_vec;
dat(~sig) = sentinel;

end
