function tfce_vec = tfce_one_fmri_dat( ...
       t_vec, removed_voxels, volInfo, ...
       voxsize, H, E, conn, sidedness, tail)
% TFCE_ONE_FMRI_DAT
% Single t-map TFCE transform (parfor-safe)

dim = volInfo.dim;
vol = zeros(dim,'double');
vol(~removed_voxels) = t_vec;
vol(~isfinite(vol)) = 0;

switch sidedness
   case 'one'
       if strcmp(tail,'pos')
           tf = call_pTFCE(max(vol,0), voxsize, H, E, conn);
       else
           tf = call_pTFCE(max(-vol,0), voxsize, H, E, conn);
       end
   case 'two'
       tf_pos = call_pTFCE(max(vol,0), voxsize, H, E, conn);
       tf_neg = call_pTFCE(max(-vol,0), voxsize, H, E, conn);
       tf = max(tf_pos, tf_neg);
end

tf(~isfinite(tf) | tf < 0) = 0;
tfce_vec = tf(~removed_voxels);

end

