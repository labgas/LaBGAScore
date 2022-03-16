%% LaBGAScore_first_s2_diagnose_firstlvl.m
%
% This script will run diagnostics on first level models and publish them
% as an html report using the CANlab function scn_spm_design_check, and
% save variance inflation factors as a .mat file
%
% DOCUMENTATION
% help scn_spm_design_check
% 
% USAGE
% Script should be called from LaBGAScore_first_s1_prep_firstlvl.m, it will
% not work independently
%
% DEPENDENCIES
% The following software needs to be on your Matlab path
% - spm12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
% - CanlabCore (https://github.com/canlab/CanlabCore)
%
%__________________________________________________________________________
%
% authors: Lukas Van Oudenhove
% date:   March, 2022
%
%__________________________________________________________________________
% @(#)% LaBGAScore_first_s2_diagnose_firstlvl.m         v1.0        
% last modified: 2022/03/14
%
%
%% RUN DIAGNOSTICS

vifs = scn_spm_design_check(subjfirstdir,'events_only','vif_thresh',3);
save('vifs','vifs');
