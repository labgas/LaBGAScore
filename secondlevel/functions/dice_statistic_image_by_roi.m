function roi_table = dice_statistic_image_by_roi(imgA, imgB, roi_img, varargin)
% DICE_STATISTIC_IMAGE_BY_ROI
%
% Compute per-ROI Dice coefficients between two statistic_image objects.
%
% INPUTS
%  imgA, imgB : statistic_image
%  roi_img   : ROI definition (statistic_image or fmri_data)
%
% OPTIONAL NAME–VALUE
%  'binary_mask' : {maskA, maskB}
%  'p_threshold' : p-value threshold (default = 0.05)
%
% OUTPUT
%  roi_table : table with columns:
%      ROI      ROI_label
%      Dice     Dice coefficient
%      nA       Voxels in A within ROI
%      nB       Voxels in B within ROI
%      nOverlap Overlapping voxels in ROI
%

%% defaults
p_thresh = 0.05;
binary_mask = [];

for i = 1:2:numel(varargin)
   switch lower(varargin{i})
       case 'p_threshold'
           p_thresh = varargin{i+1};
       case 'binary_mask'
           binary_mask = varargin{i+1};
   end
end

%% masks
if isempty(binary_mask)
   maskA = imgA.p < p_thresh;
   maskB = imgB.p < p_thresh;
else
   maskA = logical(binary_mask{1}(:));
   maskB = logical(binary_mask{2}(:));
end

%% ROI labels
roi_labels = unique(roi_img.dat);
roi_labels(roi_labels == 0 | isnan(roi_labels)) = [];

ROI = [];
Dice = [];
nA = [];
nB = [];
nOverlap = [];

for r = 1:numel(roi_labels)
   lab = roi_labels(r);
   roi_mask = roi_img.dat == lab;

   A = maskA & roi_mask;
   B = maskB & roi_mask;

   nA_r = sum(A);
   nB_r = sum(B);
   nOv = sum(A & B);

   if nA_r + nB_r == 0
       dice_r = NaN;
   else
       dice_r = 2 * nOv / (nA_r + nB_r);
   end

   ROI(end+1,1) = lab;
   Dice(end+1,1) = dice_r;
   nA(end+1,1) = nA_r;
   nB(end+1,1) = nB_r;
   nOverlap(end+1,1) = nOv;
end

roi_table = table(ROI, Dice, nA, nB, nOverlap);

end

