function [tfce_vec, dh_used] = tfce_one_fmri_dat( ...
        t_vec, removed_voxels, volInfo, H, E, conn, sidedness, tail, dh)
% tfce_one_fmri_dat  Single statistic-map TFCE transform (parfor-safe).
%
%   [tfce_vec, dh_used] = tfce_one_fmri_dat(t_vec, removed_voxels, volInfo, ...
%                                           H, E, conn, sidedness, tail, dh)
%
%   Rebuilds a 3D volume from a masked statistic vector, TFCE-transforms it,
%   and returns the values back in the masked-vector layout. Takes only plain
%   data (no objects), so it is safe to call from inside parfor.
%
%   Inputs
%     t_vec           [nInMask x 1] statistic values
%     removed_voxels  logical mask from the fmri_data object
%     volInfo         struct with .dim
%     H, E, conn      TFCE height exponent, extent exponent, connectivity
%     sidedness       'one' | 'two'
%     tail            'pos' | 'neg' (only when sidedness is 'one')
%     dh              integration step; pass [] on the observed map and reuse
%                     the returned dh_used for every permutation
%
%   Outputs
%     tfce_vec        [nInMask x 1] TFCE values
%     dh_used         the step actually used
%
%   SIGNATURE CHANGE
%     The old signature was
%         tfce_vec = tfce_one_fmri_dat(t_vec, removed_voxels, volInfo, ...
%                                      voxsize, H, E, conn, sidedness, tail)
%     `voxsize` is gone. It existed only because the previous implementation
%     forwarded it to pTFCE, which does not use TFCE's H/E parameters at all --
%     the voxel size was landing in pTFCE's `mask` slot and H/E in its resel-
%     count and voxel-count slots. Classic TFCE has no use for voxel size.
%     `dh` is new and optional.
%
%   See also tfce_volume, tfce_transform_3d, group_tfce_from_subject_maps.

if nargin < 9, dh = []; end

dim = volInfo.dim;
vol = zeros(dim,'double');
vol(~removed_voxels) = t_vec;
vol(~isfinite(vol)) = 0;

[tf, dh_used] = tfce_volume(vol, H, E, conn, sidedness, tail, dh);

tfce_vec = tf(~removed_voxels);

end
