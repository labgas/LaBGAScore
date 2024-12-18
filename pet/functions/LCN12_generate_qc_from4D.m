function LCN12_generate_qc_from4D(img4D,outputdir,name)
% 
% this function will generate 3D images in each orthogonal direction
% (axial, sagital, coronal) of the middle slice in which the third
% dimension represents the fourth dimension of the 4D image. These images
% can be used to evaluate e.g. movement if the fourth dimension represents
% time.
%
% INPUT 
%   img4D = 4D matrix
%   outputdir = optional, name of the output directory
%   name      = optional, name of the outputfile (no extension)
%
% author: Patrick Dupont
% date: June 2021
% history: 
%__________________________________________________________________________
% @(#)LCN12_generate_qc_from4D.m        v0.1      last modified: 2021/06/02

if nargin == 1
   outputdir = '.';
   name = 'qc';
elseif nargin == 2
   name = 'qc';
end

dimensions = size(img4D);

% generate images
qcimg1 = squeeze(img4D(:,:,round(dimensions(3)/2),:));
qcimg2 = squeeze(img4D(:,round(dimensions(2)/2),:,:));
qcimg3 = squeeze(img4D(round(dimensions(1)/2),:,:,:));

% generate filename
filename_qc1 = fullfile(outputdir,[name '1.nii']);
filename_qc2 = fullfile(outputdir,[name '2.nii']);
filename_qc3 = fullfile(outputdir,[name '3.nii']);

% write images to file
LCN12_write_image(qcimg1,filename_qc1,'qc image');        
LCN12_write_image(qcimg2,filename_qc2,'qc image');        
LCN12_write_image(qcimg3,filename_qc3,'qc image');        

end