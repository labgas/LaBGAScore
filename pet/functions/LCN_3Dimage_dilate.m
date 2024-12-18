function outputimage = LCN_3Dimage_dilate(inputimage,kernel)
% inputimage = 3D binary image
% kernel = vector with the number of voxels to dilate in each of the 3 dimensions
% outputimage = 3D binary image of the same size as inputimage
% 
% Needs to be extended but for what we want, this is sufficient
%__________________________________________________________________________
%
% author: Patrick Dupont
% date:   May 2018
% history: kernel size adapted
%__________________________________________________________________________
% @(#) LCN_3Dimage_dilate.m      v0.11            last modified: 2018/05/24

% make kernel in 3D
kernel3D = ones(kernel(1)+1,kernel(2)+1,kernel(3)+1);

% in case inputimage is not binary, treat all nonzero voxels as part of the
% binary image
inputimage = (inputimage ~= 0);

% convolve with kernel
outputimage = convn(inputimage,kernel3D,'same');

% binarize outputimage
outputimage = double(outputimage > 0);

end