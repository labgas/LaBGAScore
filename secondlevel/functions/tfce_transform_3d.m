function tfce = tfce_transform_3d(stat, H, E, conn, dh)
% Standard explicit TFCE (Smith & Nichols, 2009)

if nargin < 5 || isempty(dh)
   maxval = max(stat(:));
   if maxval <= 0 || isnan(maxval)
       tfce = zeros(size(stat),'single');
       return
   end
   dh = maxval / 100;  % 100 thresholds by default
end

tfce = zeros(size(stat),'single');
hvals = 0:dh:max(stat(:));
hvals = hvals(hvals > 0);

for h = hvals
   bw = stat >= h;
   CC = bwconncomp(bw, conn);
   for c = 1:CC.NumObjects
       idx = CC.PixelIdxList{c};
       tfce(idx) = tfce(idx) + (h^H) * (numel(idx)^E) * dh;
   end
end

end

