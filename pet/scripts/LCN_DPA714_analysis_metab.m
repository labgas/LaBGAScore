% LCN_DPA714_analysis_metab.m
%
% *NOTE FROM LaBGAS -- THIS IS VENDORED LEGACY CODE, NOT ADAPTED*
%
%   Kept because nothing else in this repo does what it does: it compares FOUR
%   alternative models for the parent-fraction (metabolite) correction of DPA714
%   - mono- and bi-exponential, each constrained and unconstrained - against the
%   Hill model that LaBGAScore_pet_model_TSPO_DPA714.m actually uses, and ranks
%   them by AIC and Schwarz criterion. It is how the Hill choice was arrived at.
%
%   IT IS A HAND-RUN SCRIPT IN ITS ORIGINAL KU Leuven FORM. Unlike
%   LaBGAScore_pet_model_TSPO_DPA714.m and LaBGAScore_pet_preprocess_data.m,
%   which are LaBGAS adaptations of their LCN12 ancestors, this one was never
%   adapted. Before running it you must replace, at minimum:
%
%     - the hardcoded Windows paths (lines ~21, 32-33, 37: 'C:\DATA\TSPO'),
%       which on this setup come from <study>_prep_s0_define_directories
%     - the hardcoded SUBJECTS list (line ~26)
%     - xlswrite -> writetable (lines ~313, 318-319); xlswrite needs Excel/COM
%       and does not work on Linux
%     - date -> char(datetime('today')) (lines ~316-317); date is deprecated
%
%   It calls the eight LCN_{calc,cost}_intact_tracer_{mono,bi}exp[_con|_delay]
%   functions in ../functions, which exist in this repo ONLY for this script -
%   the TSPO pipeline itself uses the _hill pair alone.
%
% -------------------------------------------------------------------------
%
% Original header follows.
%
% LCN_DPA714_analysis_metab
%
% This script will calculate different models for metabolite correction of
% DPA714 data.
%
% author: Patrick Dupont
% date: July, 2023
% history: July 2023: tested for matlab2018B and matlab2022B
%
%__________________________________________________________________________
% @(#)LCN_DPA714_analysis_metab.m       v0.1      last modified: 2023/07/11

clear
close all

global TIME_METAB
global FRACTION_INTACT_TRACER
global WEIGHTS_METAB

%------------- SETTINGS ---------------------------------------------------
maindir           = 'C:\DATA\TSPO'; % directory where the folders of each subject can  be found
sessiondir        = ''; % if empty, we assume that there is no folder session and the folders anat and pet are directly under the subject folder
infostring_tracer = 'trc-DPA714';
infostring_metab = 'data_metab'; % % in the folder pet, we assume a .m file name "subjectname"_"sessiondir"_"infostring_tracer"_"infostring_metab".m (example sub-test_trc-DPA714_rec-acdyn_data_metab.nii)

SUBJECTS  = {
% subjectdir
'sub-KUL102B'
'sub-KUL102'
};

logfile    = 'C:\DATA\TSPO\logfile_analysis_metab_data.txt';
resultsdir = 'C:\DATA\TSPO';

SUBJECT_DATA = {
    % location of metab data
    'C:\DATA\TSPO\sub-KUL102\sub-KUL102_trc-DPA714_data_metab.m'
};


figures_on = 0; % if 1, the script will pause after showing each figure and the user has to hit a key to proceed.
save_figures = 1; % figures will be saved as matlab figures. The programme will not pause in this case.
%--------------- SETTINGS METABOLITE CORRECTION ---------------------------
p0_biexp_delay = [5 0.8 40 2]; % see LCN_calc_intact_tracer_biexp_delay for details
p0_biexp_con   = [5 0.8 40];   % see LCN_calc_intact_tracer_biexp_con for details
p0_hill        = [50 -1];      % see LCN_calc_intact_tracer_hill for details
p0_monoexp     = [50 0.9];     % see LCN_calc_intact_tracer_monoexp for details
p0_monoexp_con = 50;           % see LCN_calc_intact_tracer_monoexp_con for details
% decide until which time you want to include samples used for the
% metabolite correction
% only the models in the model_list will be tested
model_list = {
    'biexp_delay'    
    'biexp_con'
    'hill'
    'monoexp'
    'monoexp_con'
};

%++++ END OF SETTINGS - DO NOT CHANGE BELOW THIS LINE +++++++++++++++++++++

nr_subjects = size(SUBJECTS,1);
nr_models = size(model_list,1);
curdir = pwd;
fid  = fopen(logfile,'a+');
fprintf(fid,'Analysis %s\n',datetime('now'));
fprintf(fid,'Settings\n');
options = optimset('Display','off','Algorithm','active-set','MaxFunEvals',1000); %'active-set', 'trust-region-reflective', 'interior-point',  'interior-point-convex', 'levenberg-marquardt', 'trust-region-dogleg',  'lm-line-search', or 'sqp'.
% initialize
outputdata_metab_AIC = cell(nr_models+1,nr_subjects+1);
outputdata_metab_SC  = cell(nr_models+1,nr_subjects+1);
outputdata_metab_AIC{1,1} = 'model/subject';
outputdata_metab_SC{1,1}  = 'model/subject';
go = 1;

for i = 1:nr_models
    if strcmp(model_list{i},'biexp_delay') == 1
       nr_params_biexp_delay = length(p0_biexp_delay);
       outputdata_metab_AIC{i+1,1} = 'biexp_delay';
       outputdata_metab_SC{i+1,1}  = 'biexp_delay';
       fprintf(fid,'MODEL BIEXP_DELAY: initial conditions (see LCN_calc_intact_tracer_biexp_delay.m) %4.2f %4.2f %4.2f %4.2f\n',p0_biexp_delay(1),p0_biexp_delay(2),p0_biexp_delay(3),p0_biexp_delay(4));
    elseif strcmp(model_list{i},'biexp_con') == 1
       nr_params_biexp_con = length(p0_biexp_con);
       outputdata_metab_AIC{i+1,1} = 'biexp_con';
       outputdata_metab_SC{i+1,1}  = 'biexp_con';
       fprintf(fid,'MODEL BIEXP_CON: initial conditions (see LCN_calc_intact_tracer_biexp_con.m) %4.2f %4.2f %4.2f \n',p0_biexp_con(1),p0_biexp_con(2),p0_biexp_con(3));
    elseif strcmp(model_list{i},'hill') == 1
       nr_params_hill = length(p0_hill);
       outputdata_metab_AIC{i+1,1} = 'hill';
       outputdata_metab_SC{i+1,1}  = 'hill';
       fprintf(fid,'MODEL HILL: initial conditions (see LCN_calc_intact_tracer_hill.m) %4.2f %4.2f \n',p0_hill(1),p0_hill(2));
    elseif strcmp(model_list{i},'monoexp') == 1
       nr_params_monoexp = length(p0_monoexp);
       outputdata_metab_AIC{i+1,1} = 'monoexp';
       outputdata_metab_SC{i+1,1}  = 'monoexp';
       fprintf(fid,'MODEL MONOEXP: initial conditions (see LCN_calc_intact_tracer_monoexp.m) %4.2f %4.2f \n',p0_monoexp(1),p0_monoexp(2));
    elseif strcmp(model_list{i},'monoexp_con') == 1
       nr_params_monoexp_con = length(p0_monoexp_con);
       outputdata_metab_AIC{i+1,1} = 'monoexp_con';
       outputdata_metab_SC{i+1,1}  = 'monoexp_con';
       fprintf(fid,'MODEL MONOEXP_CON: initial conditions (see LCN_calc_intact_tracer_monoexp_con.m) %4.2f \n',p0_monoexp_con(1));
    else
       fprintf('model %s not implementend \n',model_list{i});
       fprintf(fid,'model %s not implementend \n',model_list{i});
    end
end

for subj = 1:nr_subjects
    clear subjectdir subjectname data_metab fine_time nr_samples    
    clear outputdata_metab filename_metab hfig
    
    subjectname = SUBJECTS{subj,1};
    subjectdir  = fullfile(fullfile(maindir,subjectname),['pet_' infostring_tracer]);

    fprintf('working on subject %s \n',subjectname);

    % reset subject specific global variables
    TIME_METAB    = [];
    FRACTION_INTACT_TRACER = [];
    WEIGHTS_METAB = [];
        
    % initialize
    outputdata_metab = cell(nr_models+1,8);
    
    % find the metab file
    [filename_metab,go]   = LCN_check_filename(subjectdir,[subjectname '*_' infostring_tracer '*_*' infostring_metab '.m']);
        
    outputdata_metab_AIC{1,subj+1} = subjectname;
    outputdata_metab_SC{1,subj+1}  = subjectname;
        
    % read the file containing the variable data_metab
    %-------------------------------------------------
    fprintf(fid,'Metab datafile:  %s\n',filename_metab);
    copyfile(filename_metab,'tmp_metab_data.m');
    tmp_metab_data; % the variable frames_timing is now known
    delete('tmp_metab_data.m');

    TIME_METAB             = data_metab(:,1)/60;
    FRACTION_INTACT_TRACER = data_metab(:,2)/100;
    WEIGHTS_METAB          = data_metab(:,3);

    % normalize WEIGHTS
    WEIGHTS_METAB = WEIGHTS_METAB/sum(WEIGHTS_METAB);        
    
    fine_time  = 0:0.01:max(TIME_METAB);
    nr_samples = length(TIME_METAB);
        
    % make header for outputdata per subject
    outputdata_metab{1,1} = 'model\results';
    outputdata_metab{1,2} = 'goodness of fit';
    outputdata_metab{1,3} = 'AIC';
    outputdata_metab{1,4} = 'SC';
    outputdata_metab{1,5} = 'p1'; 
    outputdata_metab{1,6} = 'p2'; 
    outputdata_metab{1,7} = 'p3'; 
    outputdata_metab{1,8} = 'p4'; 
       
    if figures_on == 1 || save_figures == 1
       hfig = figure(subj);
       set(hfig,'Name',subjectname);
    end 

    for j = 1:nr_models
        if strcmp(model_list{j},'biexp_delay') == 1
           % biexp_delay
           %++++++++++++
           clear pfit res intact_fine
           % fit model
           pfit  = fminsearch('LCN_cost_intact_tracer_biexp_delay',p0_biexp_delay);   
           % determine the error of the fit
           [res] = LCN_cost_intact_tracer_biexp_delay(pfit);
             
           intact_fine = LCN_calc_intact_tracer_biexp_delay(pfit,fine_time);     
           if figures_on == 1 || save_figures == 1
              subplot(1,nr_models,j);
              plot(TIME_METAB,FRACTION_INTACT_TRACER,'o')
              axis([0 1.1*max(TIME_METAB) 0 1])
              hold on
              plot(fine_time,intact_fine)
              title('biexp\_delay')
              xlabel('time (min)')
              ylabel('fraction intact tracer')
           end
           outputdata_metab_AIC{j+1,subj+1} = nr_samples.*log(min(res))+2.*nr_params_biexp_delay;
           outputdata_metab_SC{j+1,subj+1}  = nr_samples.*log(min(res))+nr_params_biexp_delay.*log(nr_samples);
           outputdata_metab{1+j,1} = 'biexp_delay';
           outputdata_metab{1+j,2} = res;
           outputdata_metab{1+j,3} = outputdata_metab_AIC{j+1,subj+1};
           outputdata_metab{1+j,4} = outputdata_metab_SC{j+1,subj+1};
           outputdata_metab{1+j,5} = pfit(1);
           outputdata_metab{1+j,6} = pfit(2);
           outputdata_metab{1+j,7} = pfit(3);
           outputdata_metab{1+j,8} = pfit(4);
        elseif strcmp(model_list{j},'biexp_con') == 1
           clear pfit res intact_fine
           % biexp_con
           %++++++++++
           % fit model
           pfit  = fminsearch('LCN_cost_intact_tracer_biexp_con',p0_biexp_con);   
           % determine the error of the fit
           [res] = LCN_cost_intact_tracer_biexp_con(pfit);   
           intact_fine = LCN_calc_intact_tracer_biexp_con(pfit,fine_time);     
           if figures_on == 1 || save_figures == 1
              subplot(1,nr_models,j);
              plot(TIME_METAB,FRACTION_INTACT_TRACER,'o')
              axis([0 1.1*max(TIME_METAB) 0 1])
              hold on
              plot(fine_time,intact_fine)
              title('biexp\_con');
              xlabel('time (min)')
              ylabel('fraction intact tracer')
           end
           outputdata_metab_AIC{j+1,subj+1} = nr_samples.*log(min(res))+2.*nr_params_biexp_con;
           outputdata_metab_SC{j+1,subj+1}  = nr_samples.*log(min(res))+nr_params_biexp_con.*log(nr_samples);
           outputdata_metab{1+j,1} = 'biexp_con';
           outputdata_metab{1+j,2} = res;
           outputdata_metab{1+j,3} = outputdata_metab_AIC{j+1,subj+1};
           outputdata_metab{1+j,4} = outputdata_metab_SC{j+1,subj+1};
           outputdata_metab{1+j,5} = pfit(1);
           outputdata_metab{1+j,6} = pfit(2);
           outputdata_metab{1+j,7} = pfit(3);
        elseif strcmp(model_list{j},'hill') == 1
           clear pfit res intact_fine
           % hill
           %+++++
           % fit model
           pfit  = fminsearch('LCN_cost_intact_tracer_hill',p0_hill);   
           % determine the error of the fit
           [res] = LCN_cost_intact_tracer_hill(pfit);   
           intact_fine = LCN_calc_intact_tracer_hill(pfit,fine_time);     
           if figures_on == 1 || save_figures == 1
              subplot(1,nr_models,j);
              plot(TIME_METAB,FRACTION_INTACT_TRACER,'o')
              axis([0 1.1*max(TIME_METAB) 0 1])
              hold on
              plot(fine_time,intact_fine)
              title('hill');
              xlabel('time (min)')
              ylabel('fraction intact tracer')
           end
           outputdata_metab_AIC{j+1,subj+1} = nr_samples.*log(min(res))+2.*nr_params_hill;
           outputdata_metab_SC{j+1,subj+1}  = nr_samples.*log(min(res))+nr_params_hill.*log(nr_samples);
           outputdata_metab{1+j,1} = 'hill';
           outputdata_metab{1+j,2} = res;
           outputdata_metab{1+j,3} = outputdata_metab_AIC{j+1,subj+1};
           outputdata_metab{1+j,4} = outputdata_metab_SC{j+1,subj+1};
           outputdata_metab{1+j,5} = pfit(1);
           outputdata_metab{1+j,6} = pfit(2);
        elseif strcmp(model_list{j},'monoexp') == 1
           clear pfit res intact_fine
           % monoexp
           %+++++
           % fit model
           pfit  = fminsearch('LCN_cost_intact_tracer_monoexp',p0_monoexp);   
           % determine the error of the fit
           [res] = LCN_cost_intact_tracer_monoexp(pfit);   
           intact_fine = LCN_calc_intact_tracer_monoexp(pfit,fine_time);     
           if figures_on == 1 || save_figures == 1
              subplot(1,nr_models,j);
              plot(TIME_METAB,FRACTION_INTACT_TRACER,'o')
              axis([0 1.1*max(TIME_METAB) 0 1])
              hold on
              plot(fine_time,intact_fine)
              title('monoexp');
              xlabel('time (min)')
              ylabel('fraction intact tracer')
           end
           outputdata_metab_AIC{j+1,subj+1} = nr_samples.*log(min(res))+2.*nr_params_monoexp;
           outputdata_metab_SC{j+1,subj+1}  = nr_samples.*log(min(res))+nr_params_monoexp.*log(nr_samples);
           outputdata_metab{1+j,1} = 'monoexp';
           outputdata_metab{1+j,2} = res;
           outputdata_metab{1+j,3} = outputdata_metab_AIC{j+1,subj+1};
           outputdata_metab{1+j,4} = outputdata_metab_SC{j+1,subj+1};
           outputdata_metab{1+j,5} = pfit(1);
           outputdata_metab{1+j,6} = pfit(2);
        elseif strcmp(model_list{j},'monoexp_con') == 1
           clear pfit res intact_fine
           % monoexp_con
           %++++++++++++
           % fit model
           pfit  = fminsearch('LCN_cost_intact_tracer_monoexp_con',p0_monoexp_con);   
           % determine the error of the fit
           [res] = LCN_cost_intact_tracer_monoexp_con(pfit);   
           intact_fine = LCN_calc_intact_tracer_monoexp_con(pfit,fine_time);     
           if figures_on == 1 || save_figures == 1
              subplot(1,nr_models,j);
              plot(TIME_METAB,FRACTION_INTACT_TRACER,'o')
              axis([0 1.1*max(TIME_METAB) 0 1])
              hold on
              plot(fine_time,intact_fine)
              title('monoexp\_con');
              xlabel('time (min)')
              ylabel('fraction intact tracer')
           end
           outputdata_metab_AIC{j+1,subj+1} = nr_samples.*log(min(res))+2.*nr_params_monoexp_con;
           outputdata_metab_SC{j+1,subj+1}  = nr_samples.*log(min(res))+nr_params_monoexp_con.*log(nr_samples);
           outputdata_metab{1+j,1} = 'monoexp_con';
           outputdata_metab{1+j,2} = res;
           outputdata_metab{1+j,3} = outputdata_metab_AIC{j+1,subj+1};
           outputdata_metab{1+j,4} = outputdata_metab_SC{j+1,subj+1};
           outputdata_metab{1+j,5} = pfit(1);
        end
    end
    cd(subjectdir);
    % save results (figures and fittings)
    if figures_on == 1 
       pause;
    end
    if save_figures == 1
       saveas(hfig,[subjectname 'fig_metab.fig']);
    end
    outputfile_excel = [subjectname '_results_metab.xls'];
    xlswrite(outputfile_excel,outputdata_metab);
    eval(['save ' subjectname '_metab outputdata_metab']);     
end
outputfile_excel_AIC = fullfile(resultsdir,['Results_metab_AIC_' date '.xls']);
outputfile_excel_SC  = fullfile(resultsdir,['Results_metab_SC_' date '.xls']);
xlswrite(outputfile_excel_AIC,outputdata_metab_AIC);
xlswrite(outputfile_excel_SC,outputdata_metab_SC);
cd(resultsdir);
eval(['save Results_metab_' date '_outputdata_metab_AIC outputdata_metab_SC']);     

% close log file
%---------------
fclose(fid);
fprintf('all done \n');