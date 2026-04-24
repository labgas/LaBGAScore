function [dice, info] = dice_statistic_image(img1, img2, varargin)
%
% DICE_STATISTIC_IMAGE
%
% Compute the Dice similarity coefficient between two Canlab statistic_image objects
% by thresholding their p-values or using explicitly provided binary masks.
%
% INPUTS
%  img1, img2 : statistic_image
%
% OPTIONAL NAME–VALUE PAIRS
%  'p_threshold' : p-value threshold (default = 0.05)
%  'binary_mask' : cell array {mask1, mask2} of logical vectors
%                  (overrides p-value thresholding)
%
% OUTPUTS
%  dice : Dice similarity coefficient
%
%  info : struct with fields
%      .nA       number of suprathreshold voxels in image A
%      .nB       number of suprathreshold voxels in image B
%      .nOverlap number of overlapping suprathreshold voxels
%      .maskA    binary mask for image A
%      .maskB    binary mask for image B
%

%% ---------------- defaults ----------------
p_thresh = 0.05;
binary_mask = [];

for i = 1:2:numel(varargin)
   switch lower(varargin{i})
       case 'p_threshold'
           p_thresh = varargin{i+1};
       case 'binary_mask'
           binary_mask = varargin{i+1};
       otherwise
           error('Unknown option: %s', varargin{i});
   end
end

%% ---------------- checks ----------------
if ~isa(img1, 'statistic_image') || ~isa(img2, 'statistic_image')
   error('Inputs must be Canlab statistic_image objects.');
end

if numel(img1.dat) ~= numel(img2.dat)
   error('Images do not have matching voxel dimensions.');
end

%% :white_check_mark: Repair geometry ↔ data consistency (Canlab-native)
% try
%    img1 = trim_mask(img1);       % works if masked
%    img2 = trim_mask(img2);       % works if masked
% catch
%    % trim_mask errors if object is not maskable
% end
% 
% img1 = replace_empty(img1);       % fixes remaining small mismatches
% img2 = replace_empty(img2);       % fixes remaining small mismatches


%% ---------------- build binary masks ----------------
if ~isempty(binary_mask)
   % User-supplied masks
   if numel(binary_mask) ~= 2
       error('binary_mask must be a cell array with two elements.');
   end

   maskA = logical(binary_mask{1}(:));
   maskB = logical(binary_mask{2}(:));

else
   % Threshold using p-values
   if isempty(img1.p) || isempty(img2.p)
       error('Images do not contain p-values. Provide binary_mask explicitly.');
   end

   % valid voxels: not NaN and not masked out
   valid = ~(isnan(img1.p) | isnan(img2.p));

   maskA = false(size(img1.p));
   maskB = false(size(img2.p));

   maskA(valid) = img1.p(valid) <= p_thresh;
   maskB(valid) = img2.p(valid) <= p_thresh;
end

%% ---------------- Dice coefficient ----------------
nA = sum(maskA);
nB = sum(maskB);
nOverlap = sum(maskA & maskB);

if nA + nB == 0
   dice = NaN; % undefined if neither map has suprathreshold voxels
else
   dice = 2 * nOverlap / (nA + nB);
end

%% ---------------- output info ----------------
info = struct();
info.nA = nA;
info.nB = nB;
info.nOverlap = nOverlap;
info.maskA = maskA;
info.maskB = maskB;

end

