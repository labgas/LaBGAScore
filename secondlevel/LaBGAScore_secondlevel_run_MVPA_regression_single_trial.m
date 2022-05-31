%%% LaBGAScore_secondlevel_run_MVPA_regression_single_trial
%
% This script runs MVPA regression analysis on a continuous outcome Y
% (lasso-pcr, but can easily be adapted to support vector regression or
% other methods available in CANlab's predict function) on an fmri_data_st
% object created using LaBGAS_secondlevel_create_single_trial_fmri_data_st_obj.m. 
% That script should be run first, or the present script will load the data
% object if it is saved by the previous script.
%
% It is based on the extremely helpful tutorials on single trial analysis
% in the context of MVPA by @bogpetre @CANlab.
%
% Here are Bogdan's walkthroughs:
% https://canlab.github.io/_pages/canlab_single_trials_demo/demo_norming_comparison.html
% https://canlab.github.io/_pages/mlpcr_demo/mlpcr_demo.html (WiP)
%
% Here are two scripts @lukasvo76 adapted from these walkthroughs
% https://www.dropbox.com/sh/e17nl3ew1db1twk/AACO9QAEt6Sy3TejH-n-tbdEa?dl=0
% https://www.dropbox.com/sh/bm0at2dr81isk70/AABD67D_bF8A0NFa4gtt2dHNa?dl=0
% 
% Another highly helpful resource in this context is this Nature Methods
% paper by Tor and Wani Woo
% https://www.nature.com/articles/s41596-019-0289-5
%
% @lukasvo76's version of the script for this paper can be found here
% https://www.dropbox.com/sh/v2nsgoqmbi0cqnk/AAD6I1Gn5KUM6aViom4TLeVJa?dl=0
%
% NOTE: this script may be replaced/complemented by a new one based on
% Bogdan's more recent object-oriented machine learning framework for fMRI
% data
% 
% NOTE: this script is work in progress, not yet extensively tested
% 
%__________________________________________________________________________
%
% author: lukas.vanoudenhove@kuleuven.be
% date:   April, 2021
%__________________________________________________________________________
% @(#)% LaBGAScore_secondlevel_run_MVPA_regression_single_trial     v2.0        
% last modified: 2021/05/31


%% LOAD FMRI_DATA_ST OBJECT AND OTHER NECESSARY VARIABLES IF NEEDED
%--------------------------------------------------------------------------

if ~exist('resultsdir','var')
    a_set_up_paths_always_run_first
end

if ~exist('DSGN','var') || ~exist('DAT','var')
    load(fullfile(resultsdir,'image_names_and_setup.mat'));
    if ~isfield(DAT,'BEHAVIOR')
        error('\n Behavioral data not yet added to DAT structure - run prep_1b script first')
    end
end

if ~exist('fmri_dat','var')
    load(fullfile(resultsdir,['single_trial_fmri_data_st_object_' DSGN.modelingfilesdir '.mat']));
end


%% DEFINE SUBJECT IDENTIFIERS
%--------------------------------------------------------------------------

subject_id = fmri_dat.metadata_table.(subj_identifier);
[uniq_subject_id, ~, subject_id] = unique(subject_id,'stable');
n_subj = size(uniq_subject_id,1);

%% MASK AND Z-SCORE DV IF REQUESTED IN OPTIONS
%--------------------------------------------------------------------------

% MASK

if ~isempty(maskname_mvpa_reg_st)
    mask = fmri_mask_image(maskname_mvpa_reg_st);
    fmri_dat = fmri_dat.apply_mask(mask);
    fmri_dat.mask = mask; % fmri_data.apply_mask does not seem to update mask info of the object automatically, so we do that manually here
    fmri_dat.mask_descrip = maskname_mvpa_reg_st;
end

% ZSCORE BEHAVIORAL OUTCOME
% NOTE: useful for more interpretable values of prediction MSE

if zscore_outcome
    fmri_dat.Y = zscore(fmri_dat.Y);
end


%% DATA VISUALISATION PRIOR TO MODEL BUILDING
%--------------------------------------------------------------------------

% BETA IMAGES
h1=figure;
for sub = 1:n_subj
    subj_idx = sub == subject_id;
    this_subj_dat = fmri_dat.dat(:,subj_idx);
    q(sub,:) = quantile(this_subj_dat(:),[0.025,0.5,0.975]);
    mu = mean(mean(this_subj_dat(:)));
    sd = std(this_subj_dat(:));
    h1 = plot([mu-sd, mu+sd],[sub,sub],'-');
    hold on;
    h2 = plot(mu,sub,'o');
    h2.Color = h1.Color;
end
box off
title('Distribution of beta weights');
xlabel('\beta');
ylabel('Subject');
hold off

p = get(gcf,'Position');
set(gcf,'Position',[p(1:2),1024,2048],'WindowState','Maximized');

clear sub

% BEHAVIORAL OUTCOME
% over subjects
b1=figure;
hold off;
b1=histogram(fmri_dat.Y);
box off
title(['Histogram of single trial ' behav_outcome]);
xlabel(behav_outcome);
ylabel('n(observations)');
set(gcf,'WindowState','Maximized');

% per subject
b2=figure;
for sub = 1:n_subj
    this_idx_Y = find(sub == subject_id);
    this_Y = fmri_dat.Y(this_idx_Y);

    subplot(ceil(sqrt(n_subj)), ceil(n_subj/ceil(sqrt(n_subj))), sub);
    hold off
    b2 = histogram(this_Y);
    box off
    title(uniq_subject_id{sub});
    xlabel(behav_outcome);
    ylabel('n(obs)');
end
set(gcf,'WindowState','Maximized');

clear sub


%% CROSS-VALIDATION FOLD SELECTION
%--------------------------------------------------------------------------
% NOTE: balancing over groups, stratifying over subjects (i.e. leave whole
% subject out)

if ~isempty(DAT.BETWEENPERSON.group)
    group = fmri_dat.metadata_table.(group_identifier);
    cv = cvpartition2(group, 'Group',subject_id, 'GroupKFold', 3);
        fold_labels = zeros(size(fmri_dat.dat,2),1);
        for sub = 1:cv.NumTestSets
            fold_labels(cv.test(sub)) = sub;
        end
        
else
    cv = cvpartition2(size(fmri_dat.dat,2),'Group',subject_id, 'GroupKFold', 3);
        fold_labels = zeros(size(fmri_dat.dat,2),1);
        for sub = 1:cv.NumTestSets
            fold_labels(cv.test(sub)) = sub;
        end
        
end
    
    
%% FIT SINGLE-LEVEL MVPA MODELS ON RAW AND TRANSFORMED DATA
%--------------------------------------------------------------------------
% NOTE: we use support vector regression here, but can easily be changed to
% another algorithm for continuous outcomes included in CANlab's predict
% function, for example LASSO-PCR or PLS
% type help predict in Matlab command window for more info

    switch myscaling_mvpa_reg_st

        case 'raw'

            % default
            t0 = tic;
            [d_cverr, d_stats, d_optout] = predict(fmri_dat, 'algorithm_name', 'cv_svr', ...
                'nfolds', fold_labels, 'error_type', 'mse', 'useparallel', 0, 'verbose', 0);
            d_t = toc(t0);

            fprintf('PCR r = %0.3f\n', corr(d_stats.yfit, fmri_dat.Y));

            figure
            line_plot_multisubject(fmri_dat.Y, d_stats.yfit, 'subjid', subject_id);
            xlabel({['Observed ' behav_outcome],'(average over conditions)'}); ylabel({['SVR Estimated ' behav_outcome],'(cross validated)'})

            figure
            d_stats.weight_obj.montage;

        case 'centerimages'

            fmri_dat = fmri_dat.rescale('centerimages');    

            % centered
            t0 = tic;
            [c_cverr, c_stats, c_optout] = predict(fmri_dat, 'algorithm_name', 'cv_svr', ...
                'nfolds', fold_labels, 'error_type', 'mse', 'useparallel', 0, 'verbose', 0);
            c_t = toc(t0);

            fprintf('PCR r = %0.3f\n', corr(c_stats.yfit, fmri_dat.Y));

            figure
            line_plot_multisubject(fmri_dat.Y, c_stats.yfit, 'subjid', subject_id);
            xlabel({['Observed ' behav_outcome],'(average over conditions)'}); ylabel({['SVR Estimated ' behav_outcome],'(cross validated)'})

            figure
            c_stats.weight_obj.montage;

        case 'l2norm_images'

            fmri_dat = fmri_dat.rescale('l2norm_images');

            % l2normed
            t0 = tic;
            [l2_cverr, l2_stats, l2_optout] = predict(fmri_dat, 'algorithm_name', 'cv_svr', ...
                'nfolds', fold_labels, 'error_type', 'mse', 'useparallel', 0, 'verbose', 0);
            l2_t = toc(t0);

            fprintf('PCR r = %0.3f\n', corr(l2_stats.yfit, fmri_dat.Y));

            figure
            line_plot_multisubject(fmri_dat.Y, l2_stats.yfit, 'subjid', subject_id);
            xlabel({['Observed ' behav_outcome],'(average over conditions)'}); ylabel({['SVR Estimated ' behav_outcome],'(cross validated)'})

            figure
            l2_stats.weight_obj.montage;

        case 'zscore_images'

            fmri_dat_z = fmri_dat.rescale('zscoreimages');

            % zscored
            t0 = tic;
            [z_cverr, z_stats, z_optout] = predict(fmri_dat, 'algorithm_name', 'cv_svr', ...
                'nfolds', fold_labels, 'error_type', 'mse', 'useparallel', 0, 'verbose', 0);
            z_t = toc(t0);

            fprintf('PCR r = %0.3f\n', corr(z_stats.yfit, fmri_dat.Y));

            figure
            line_plot_multisubject(fmri_dat.Y, z_stats.yfit, 'subjid', subject_id);
            xlabel({['Observed ' behav_outcome],'(average over conditions)'}); ylabel({['SVR Estimated ' behav_outcome],'(cross validated)'})

            figure
            z_stats.weight_obj.montage;

        otherwise 
            error('\ninvalid scaling option %s specified in myscaling_mvpa_reg_st variable defined in a2_set_default_options CANlabhelpexamples script, please correct\n', myscaling_mvpa_reg_st)

    end % switch scaling


%% FIT MULTILEVEL MVPA MODEL
%--------------------------------------------------------------------------

% DETERMINE MAXIMUM NUMBER OF COMPONENTS
max_comp = floor(size(fmri_dat.dat,2).*0.75 - n_subj);
% NOTE: we specify the maximum number of components as < the number of columns in
% fmri_dat.dat (n_subjects*n_conditions(in every fold)) to avoid overfitting in multilevel models, 
% where we need to leave df for the random intercepts (upper bound 1 df per random intercept hence subject)

% FIT SINGLE LEVEL MODEL
[pcr_cverr, pcr_stats,pcr_optout] = fmri_dat.predict('algorithm_name','cv_pcr',...
    'nfolds',fold_labels, 'numcomponents',max_comp);
fprintf('PCR r = %0.3f\n', corr(pcr_stats.yfit,fmri_dat.Y));

figure
line_plot_multisubject(fmri_dat.Y, pcr_stats.yfit, 'subjid', subject_id);
xlabel({['Observed ' behav_outcome],'(average over conditions)'}); ylabel({['PCR Estimated ' behav_outcome],'(cross validated)'})

figure
pcr_stats.weight_obj.montage;

% FIT MULTILEVEL MODEL W/FIXED PARAMETER ESTIMATION
% split maximum amount of components in between and within
n_bt_comp = floor(0.75*n_subj);
n_wi_comp = max_comp - n_bt_comp;
% NOTE: max between = n_subj IN EVERY FOLD (hence n_subj - 20% in 5-fold CV), 
% and you want to put more money on within since this typically explains
% more variance

% overall model prediction
[mlpcr_cverr, mlpcr_stats, mlpcr_optout] = fmri_dat.predict('algorithm_name','cv_mlpcr',...
    'nfolds',fold_labels,'numcomponents',[n_bt_comp, n_wi_comp],'subjIDs',subject_id);
fprintf('multilevel PCR r = %0.3f\n',corr(mlpcr_stats.yfit, fmri_dat.Y));
% lukasvo76: algorithm option 'cv_mlpcr' requires subject identifier,
% which makes sense since this is a multilevel/mixed model
% note that fold labels are the same, since they respect subject membership

figure
line_plot_multisubject(fmri_dat.Y, mlpcr_stats.yfit, 'subjid', subject_id);
xlabel({['Observed ' behav_outcome],'(average over conditions)'}); ylabel({['MLPCR Estimated ' behav_outcome],'(cross validated)'})

figure
mlpcr_stats.weight_obj.montage;

figure
subplot(1,2,1)
line_plot_multisubject(pcr_stats.yfit, mlpcr_stats.yfit, 'subjid', subject_id);
xlabel({'PCR model prediction'}); ylabel('Multilevel PCR model prediction');
axis square
subplot(1,2,2);
plot(pcr_optout{1}(:),mlpcr_optout{1}(:),'.');
lsline;
xlabel('PCR model weights'); ylabel('Multilevel PCR model weights');
axis square
% NOTE: contrary to @bogpetre's walkthrough, the
% pcr and the multilevel pcr models are not exactly equivalent anymore
% since I have been specifying the number of components

% get the variance explained by the between and within component
% NOTE: These functions call the same thing under the hood, 
% but simply perform cross validation using ONLY between or within
% subject models.
[mlpcr_bt_cverr, mlpcr_bt_stats] = fmri_dat.predict('algorithm_name','cv_mlpcr_bt',...
    'nfolds',fold_labels,'numcomponents',[n_bt_comp, n_wi_comp],'subjIDs',subject_id, 'verbose', 1, 'useparallel', 1);
pred_bt = mlpcr_bt_stats.yfit;

[mlpcr_wi_cverr, mlpcr_wi_stats] = fmri_dat.predict('algorithm_name','cv_mlpcr_wi',...
    'nfolds',fold_labels,'numcomponents',[n_bt_comp, n_wi_comp],'subjIDs',subject_id, 'verbose', 1, 'useparallel', 1);
pred_wi = mlpcr_wi_stats.yfit;
% lukasvo76: algorithm options are created by bogpetre

fprintf('Between subject PCR components r = %0.3f\n', corr(mlpcr_bt_stats.yfit, fmri_dat.Y));
fprintf('Within subject PCR components r = %0.3f\n', corr(mlpcr_wi_stats.yfit, fmri_dat.Y));

figure
subplot(1,2,1)
line_plot_multisubject(fmri_dat.Y, pred_bt, 'subjid', subject_id);
xlabel({['Observed ' behav_outcome]}); ylabel('Between subject components'' prediction');
axis square
subplot(1,2,2)
line_plot_multisubject(fmri_dat.Y, pred_wi, 'subjid', subject_id);
xlabel({['Observed ' behav_outcome]}); ylabel('Within subject components'' prediction');
axis square

% FIT MULITLEVEL MODEL W/ MIXED PARAMETER ESTIMATION
% NOTE: this is not part of the walkthrough yet, but @bogpetre pushed
% mlpcr3 function to CanlabCore
% bogpetre:
% main function is mlpcr3, so help mlpcr3 for usage options.
% basically it's the same as cv_mlpcr, except there's a randInt, randSlope and fitlmeOpts now
% fitlmeOpts get passed on to fitlme. It picks some sensible defaults
% randSlope makes things much slower. randInt is roughly the sae order of magnitude as running the fixed effects version
% the function defaults to fixed effects by default, so it's a drop in replacement
% for mlpcr2.m (aka cv_mlpcr)

% overall model prediction including random intercept only
[mlpcr3_cverr, mlpcr3_stats, mlpcr3_optout] = fmri_dat.predict('algorithm_name','cv_mlpcr3',...
    'nfolds',fold_labels,'numcomponents',[n_bt_comp, n_wi_comp],'subjIDs',subject_id, 'randInt', 1);
fprintf('multilevel PCR r = %0.3f\n',corr(mlpcr3_stats.yfit, fmri_dat.Y));
% NOTE: compare with code in previous section and note change of algorithm_name option to cv_mlpcr3 and addition of randInt
% option - see help mlpcr3 for more details

figure
line_plot_multisubject(fmri_dat.Y, mlpcr3_stats.yfit, 'subjid', subject_id);
xlabel({['Observed ' behav_outcome],'(average over conditions)'}); ylabel({['MLPCR3 Estimated ' behav_outcome],'(cross validated)'})

figure
mlpcr3_stats.weight_obj.montage;

figure
subplot(1,2,1)
line_plot_multisubject(mlpcr_stats.yfit, mlpcr3_stats.yfit, 'subjid', subject_id);
xlabel({'Multilevel PCR model prediction'}); ylabel('Multilevel PCR model random int model prediction');
axis square
subplot(1,2,2);
plot(mlpcr_optout{1}(:),mlpcr3_optout{1}(:),'.');
lsline;
xlabel('Multilevel PCR model weights'); ylabel('Multilevel PCR model random int model weights');
axis square
% NOTE: note that models are very similar but not exactly equivalent

% overall model prediction including random intercept and random slope
[mlpcr3rs_cverr, mlpcr3rs_stats, mlpcr3rs_optout] = fmri_dat.predict('algorithm_name','cv_mlpcr3',...
    'nfolds',fold_labels,'numcomponents',[n_bt_comp, n_wi_comp],'subjIDs',subject_id, 'randInt', 1, 'randSlope', 1);
fprintf('multilevel PCR r = %0.3f\n',corr(mlpcr3rs_stats.yfit, fmri_dat.Y));

figure
line_plot_multisubject(fmri_dat.Y, mlpcr3rs_stats.yfit, 'subjid', subject_id);
xlabel({['Observed ' behav_outcome],'(average over conditions)'}); ylabel({['MLPCR3 Estimated ' behav_outcome],'(cross validated)'})

figure
mlpcr3rs_stats.weight_obj.montage;


%% SAVE STATS FOR ALL MODELS
%--------------------------------------------------------------------------

savefilename = fullfile(resultsdir, 'single_trial_MVPA_results.mat');
save(savefilename, 'cv','d_stats','pcr_stats','mlpcr_stats','mlpcr3_stats','mlpcr3rs_stats', '-v7.3');

