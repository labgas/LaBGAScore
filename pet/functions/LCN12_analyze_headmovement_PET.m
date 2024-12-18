function bad_frames = LCN12_analyze_headmovement_PET(rp_file,threshold_translation,threshold_rotation)
% bad_frames = LCN12_analyze_headmovement_PET(rp_file,threshold_translation,threshold_rotation)
%
% This routine displays overall displacements as 
% function of scan (time).
% A log file is created with all results.
% This routine expects the realigment parameter file of SPM12 as input
% 
% INPUT 
%   rp_file = the name of the realignment file
%   threshold_translation = maximum overall movement in mm
%   threshold_rotation    = maximum overall rotation in degrees
% 
% OUTPUT
%   bad_frames = array of the frames which are exceeding the thresholds
%
% Important:
%  - requires the installation of SPM12 - http://www.fil.ion.ucl.ac.uk/spm/
%  - the path of SPM and this routine should be included in the Matlab path
%
% The package contains software (SPM12) developed under the auspices of The
% Wellcome Department of Imaging Neuroscience, a department of the
% Institute of Neurology at University College London. The copyright of
% this software remains with that of SPM12, see
% http://www.fil.ion.ucl.ac.uk/spm/.
%    
% This routine is supplied as is. 
% Comments or questions can be send to:
% Patrick.Dupont@med.kuleuven.be
%
% IMPORTANT REMARKS: 
%   - this is research software.
%
% authors: Patrick Dupont
% date: February 2019
% history: July 2019: added additional input elements
%          December 2019: minimal input parameteres and output of bad
%                         frames.
%          May 2021: small bug fix for bad frames
%
%__________________________________________________________________________
% @(#)LCN12_analyze_headmovement_PET.m     v0.21  last modified: 2021/05/28

data = load(rp_file);
nr_data_points = size(data,1);
[pth,name,~]   = fileparts(rp_file);
outputnamefig  = [fullfile(pth,name) '.fig']; 
outputnamelog  = [fullfile(pth,name) '.log']; 
    
% replace _ by \_ to make the title correct
name2 = strrep(name,'_','\_');    

% save input and results as log file
Title         = 'LCN12_analyze_headmovement_PET.m';
fid = fopen(outputnamelog,'a');
fprintf(fid,'%c','='*ones(1,50));
fprintf(fid,'\n');
fprintf(fid,'m-file:           %s\n',Title);
fprintf(fid,'Date:             %s\n',date);
fprintf(fid,'Time:             %s\n',datestr(now,13));
fprintf(fid,'realignment parameter file: %s\n',rp_file);
fprintf(fid,'threshold for translation:  %3.2f mm\n',threshold_translation);
fprintf(fid,'threshold for rotation:     %3.2f degrees\n',threshold_rotation);
hfig = figure;
data_tx    = data(:,1); 
data_ty    = data(:,2); 
data_tz    = data(:,3);
data_roll  = (data(:,4))*180/pi; % in degrees 
data_pitch = (data(:,5))*180/pi; 
data_yaw   = (data(:,6))*180/pi;
           
subplot(3,2,1)
plot(data_tx,'.:')
hold on
index_bad1 = find(abs(data_tx) > threshold_translation);
plot([1 nr_data_points],[threshold_translation threshold_translation],'r')
plot([1 nr_data_points],[-threshold_translation -threshold_translation],'r')
ylabel('X translation in mm')
title(name2)
       
subplot(3,2,3)
plot(data_ty,'.:')
hold on
index_bad2 = find(abs(data_ty) > threshold_translation);
plot([1 nr_data_points],[threshold_translation threshold_translation],'r')
plot([1 nr_data_points],[-threshold_translation -threshold_translation],'r')
ylabel('Y translation in mm')

subplot(3,2,5)
plot(data_tz,'.:')
hold on
index_bad3 = find(abs(data_tz) > threshold_translation);
plot([1 nr_data_points],[threshold_translation threshold_translation],'r')
plot([1 nr_data_points],[-threshold_translation -threshold_translation],'r')
ylabel('Z translation in mm')
       
subplot(3,2,2)
plot(data_roll,'.:')
hold on
index_bad4 = find(abs(data_roll) > threshold_rotation);
plot([1 nr_data_points],[threshold_rotation threshold_rotation],'r')
plot([1 nr_data_points],[-threshold_rotation -threshold_rotation],'r')
ylabel('roll in degrees')

subplot(3,2,4)
plot(data_pitch,'.:')
hold on
index_bad5 = find(abs(data_pitch) > threshold_rotation);
plot([1 nr_data_points],[threshold_rotation threshold_rotation],'r')
plot([1 nr_data_points],[-threshold_rotation -threshold_rotation],'r')
ylabel('pitch in degrees')

subplot(3,2,6)
plot(data_yaw,'.:')
hold on
index_bad6 = find(abs(data_yaw) > threshold_rotation);
plot([1 nr_data_points],[threshold_rotation threshold_rotation],'r')
plot([1 nr_data_points],[-threshold_rotation -threshold_rotation],'r')
xlabel('scan number')
ylabel('yaw in degrees')

bad_frames = union(union(union(union(union(index_bad1,index_bad2),index_bad3),index_bad4),index_bad5),index_bad6);

% save figure
saveas(hfig,outputnamefig);
fclose(fid);

return;        