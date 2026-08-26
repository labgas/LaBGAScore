function tfce = tfce_transform_3d(stat, H, E, conn, dh)
% tfce_transform_3d  Threshold-free cluster enhancement (Smith & Nichols 2009).
%
%   tfce = tfce_transform_3d(stat, H, E, conn, dh)
%
%   Computes, for each voxel v,
%
%       TFCE(v) = sum over thresholds h of  h^H * extent(h,v)^E * dh
%
%   where extent(h,v) is the size of the supra-threshold cluster containing v
%   at height h, found with bwconncomp.
%
%   Inputs
%     stat  3D statistic volume. Only POSITIVE thresholds are integrated (see
%           below), so pass max(vol,0) or max(-vol,0) for the tail you want --
%           tfce_volume does this for you.
%     H     height exponent (2 is the standard default)
%     E     extent exponent (0.5 is the standard default)
%     conn  connectivity for bwconncomp (6 | 18 | 26)
%     dh    integration step. Omit or pass [] to use max(stat)/100.
%
%   Output
%     tfce  single-precision volume, same size as stat
%
%   Sign convention (previously undocumented)
%     The integration runs over h > 0 only, so negative values contribute
%     nothing. This is deliberate: a signed map is handled by transforming each
%     tail separately. Callers that pass a raw signed map get the positive tail
%     only, silently.
%
%   Sharing dh across permutations
%     The default dh depends on the map's own maximum, so two maps are
%     integrated on different grids. That is fine for a single map, but inside a
%     permutation loop it makes the null slightly incomparable to the observed
%     statistic. Derive dh once from the observed map and pass it to every
%     permutation; tfce_volume does this.
%
%   Requires the Image Processing Toolbox (bwconncomp).
%
%   See also tfce_volume, tfce_one_fmri_dat, bwconncomp.

maxval = max(stat(:));

if ~isfinite(maxval) || maxval <= 0
    tfce = zeros(size(stat),'single');
    return
end

if nargin < 5 || isempty(dh)
    dh = maxval / 100;  % 100 thresholds by default
end

if ~isfinite(dh) || dh <= 0
    error('tfce_transform_3d:badStep', 'dh must be a positive finite scalar.');
end

tfce = zeros(size(stat),'single');
hvals = dh:dh:maxval;

for h = hvals
   bw = stat >= h;
   CC = bwconncomp(bw, conn);
   for c = 1:CC.NumObjects
       idx = CC.PixelIdxList{c};
       tfce(idx) = tfce(idx) + (h^H) * (numel(idx)^E) * dh;
   end
end

end
