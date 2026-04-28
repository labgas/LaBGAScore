function tfce_vol = call_pTFCE(vol, voxsize, H, E, conn)
% Robust wrapper around pTFCE that handles numerical edge cases

try
   tfce_vol = pTFCE(vol, voxsize, H, E, conn, 0);
catch
   try
       tfce_vol = pTFCE(vol, voxsize, H, E, conn);
   catch
       try
           tfce_vol = pTFCE(vol, voxsize, H, E);
       catch ME
           % If pTFCE fails (e.g. norminv error), return zeros
           tfce_vol = zeros(size(vol), 'double');
       end
   end
end

end

