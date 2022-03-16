%% LaBGAScore_firstlevel_s3_diagnose_model.m
%
% This script will run diagnostics on first level models and publish them
% as an html report using the CANlab function scn_spm_design_check, and
% save variance inflation factors as a .mat file
%
% DOCUMENTATION
% help scn_spm_design_check
% 
% USAGE
% Script should be called from LaBGAScore_firstlevel_s1_prep_firstlvl.m, it
% is not for standalone use
%
%__________________________________________________________________________
%
% authors: Lukas Van Oudenhove
% date:   March, 2022
%
%__________________________________________________________________________
% @(#)% LaBGAScore_firstlevel_s3_diagnose_model.m         v1.0        
% last modified: 2022/03/16


%% RUN DIAGNOSTICS

vifs = scn_spm_design_check(subjfirstdir,'events_only','vif_thresh',vif_thresh);
save('vifs','vifs');