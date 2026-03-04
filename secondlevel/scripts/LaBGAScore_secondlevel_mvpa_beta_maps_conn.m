%% LaBGAScore_secondlevel_mvpa_beta_maps_conn.m
%
%
% *USAGE*
%
% This script runs MVPA regression analysis on a continuous outcome Y
% (default pcr, but can easily be adapted to pls or
% other machine learning algorithms) on an fmri_data
% object created from seed-based resting-state connectivity maps generated
% by the CONN toolbox
%
% Run this script with Matlab's publish function to generate html report of results:
% publish('LaBGAScore_mvpa_beta_maps_conn','outputDir',htmlsavedir)
%
%
% *NOTE*
% 
% This script is not yet adapted to standard dataset organization on the
% LaBGAS server, but this can easily be adapted.
%
% 
% *TUTORIALS AND DOCUMENTATION*
%
% A. CANLAB'S PREDICT FUNCTION
%
% This is the classic CANlab method of running ML models, and can be chosen
% by setting the ml_method_mvpa_reg_st option in a2_set_default_options.m
% to 'predict'
%
% This script is based on the extremely helpful tutorials on 
% single trial MVPA analysis by @bogpetre @CANlab.
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
% B. BOGDAN'S MACHINE LEARNING TOOLKIT FOR FMRI_DATA OBJECTS
%
% This is a newer method inspired by Python's scikit-learn, including more
% flexible options for algorithm and feature selection, 
% hyperparameter optimization, nested cross-validation, etc. However, it
% does require more advanced programming skills and understanding the logic
% of the method
%
% Dependency: https://github.com/canlab/ooFmriDataObjML
%
% Tutorial: https://canlab.github.io/_pages/canlab_pipelines_walkthrough/estimateBestRegionPerformance.html
% Example script: https://github.com/labgas/LaBGAScore/blob/main/secondlevel/LaBGAScore_secondlevel_ooFmriDataObjML_example.m
% 
%
% *OPTIONS*
%
%       * ml_method                     'predict' or 'oofmridataobj'
%
%           1. 'oofmridataobj': use @bogpetre's object-oriented method
%           <https://github.com/canlab/ooFmriDataObjML>
%           
%           2. 'predict': use CANlab's predict function
%           <https://github.com/canlab/CanlabCore/blob/master/CanlabCore/%40fmri_data/predict.m>
%
%       * holdout_set_method            'group' or 'no_group'
%
%           1. group: use behav_group_varname variable to balance holdout sets over groups                         
%
%           2. no_group: no group factor, stratifies by subject (i.e.leave whole subject out) since data is purely between-subject
%
%       * nfolds                        number of cross-validation folds for kfold
%
%       * algorithm_name                e.g. 'cv_pcr', or other option passed into predict function (help fmri_data.predict for options)
%
% For more extensive information about these options, see also the
% documentation for the following scripts in the CANlab_help_examples
% Github repo (LaBGAS fork)
%   <https://github.com/labgas/CANlab_help_examples/blob/master/Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/prep_3a_run_second_level_regression_and_save.m>
%   <https://github.com/labgas/CANlab_help_examples/blob/master/Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/c2f_run_MVPA_regression_single_trial.m>
%
%
% -------------------------------------------------------------------------
%
% written by: Lukas Van Oudenhove
%
% date:   KU Leuven, September, 2023
%
% -------------------------------------------------------------------------
%
% LaBGAScore_secondlevel_mvpa_beta_maps_conn.m              v1.0
%
% last modified: 2023/10/25
%
%
%
%
%% SET OPTIONS
% -------------------------------------------------------------------------

% MACHINE LEARNING OPTIONS
ml_method = 'oofmridataobj';
holdout_set_method = 'group';
nfolds = 5;
algorithm_name = 'cv_pcr';

% MASK OPTIONS
maskname = which('tpl-MNI152NLin2009cAsym_res-01_label-GM_probseg_0_15.nii.gz');


%% DEFINE PATHS, BEHAVIORAL VARIABLES, AND SEEDS
% -------------------------------------------------------------------------

% PATHS
rootdir = pwd;
inputdir = fullfile(rootdir, 'input_machine_learning');
resultsdir = fullfile(rootdir, 'results');
    if ~exist(resultsdir,'dir')
       mkdir(resultsdir);
    end
htmlsavedir = fullfile(resultsdir, 'html');
    if ~exist(htmlsavedir,'dir')
       mkdir(htmlsavedir);
    end
    
% BEHAVIORAL DATA
    
behav_filename = 'rsfMRI_IWT_conn_therapy_7_2023.xlsx';
behav_outcome_varname = 'POST_minus_PRE_NRS_ALLES';
behav_qc_varname = 'treatment_SPSS_without_exclusions';
behav_group_varname = 'treatment_SPSS';

% SEED INFO

nr_seeds = 6;
source_indices_conn = [134:139];
seednames = {'hippocampus_R','hippocampus_L','amygdala_R','amygdala_L','accumbens_R','accumbens_L'};

% MASK PREP

[~,maskname_short] = fileparts(maskname);
    if contains(maskname_short,'nii')
        [~,maskname_short] = fileparts(maskname_short);
    end
mask_string = sprintf('masked with %s', maskname_short);
mask_obj = fmri_mask_image(maskname, 'noverbose'); 
% mask_obj.dat(mask_obj.dat >= 0.15) = 1; % binarize mask
% mask_obj.dat(mask_obj.dat < 0.15) = 0; % binarize mask


%% DEFINE HELPER FUNCTION CALLED LATER
% -------------------------------------------------------------------------

dashes = '----------------------------------------------';
printstr = @(dashes) disp(dashes);
printhdr = @(str) fprintf('%s\n%s\n%s\n', dashes, str, dashes);


%% LOAD BETA MAPS AND BEHAVIORAL DATA FILE
% -------------------------------------------------------------------------

% BETAS
betas = dir(fullfile(inputdir,'BETA*'));
beta_names = {betas(:).name}';
betas_per_seed = cell(1,nr_seeds);
    for seed = 1:nr_seeds
        betas_per_seed{seed} = beta_names(contains(beta_names,num2str(source_indices_conn(seed))));
    end

% BEHAVIORAL DATA FILE    
behav_dat = readtable(fullfile(inputdir,behav_filename));

% GET RID OF SUBJECTS WITH NAN ON OUTCOME VAR
idx_nan = ~isnan(behav_dat.(behav_outcome_varname));
behav_dat = behav_dat(idx_nan,:);
    for seed = 1:nr_seeds
        betas_per_seed{seed} = betas_per_seed{seed}(idx_nan,:);
    end
    
% GET RID OF SUBJECTS EXCLUDED BASED ON qc
idx_qc = ~isnan(behav_dat.(behav_qc_varname));
behav_dat = behav_dat(idx_qc,:);
    for seed = 1:nr_seeds
        betas_per_seed{seed} = betas_per_seed{seed}(idx_qc,:);
    end
    

%% CREATE FMRI_DATA_ST OBJECTS
% -------------------------------------------------------------------------

fmri_data_objs = cell(1,nr_seeds);
Y = behav_dat.(behav_outcome_varname);
    parfor seed = 1:nr_seeds
        fmri_data_objs{seed} = fmri_data(fullfile(inputdir,betas_per_seed{seed}),mask_obj);
        fmri_data_objs{seed}.mask_descrip = maskname_short;
        fmri_data_objs{seed}.Y = Y;
        fmri_data_objs{seed}.Y_names = behav_outcome_varname;
        fmri_data_objs{seed}.dat_descrip = ['functional connectivity with seed ' seednames{seed}];
        fmri_data_objs{seed}.metadata_table = behav_dat;
        fmri_data_objs{seed}.source_notes = 'beta images from seed-based connectivity analysis in CONN';
    end
    
    
%% RUN MVPA REGRESSION MODELS
% -------------------------------------------------------------------------

mvpa_stats_results = cell(1,nr_seeds);
mvpa_stats_weight_objs = cell(1,nr_seeds);
mvpa_dats = cell(1,nr_seeds);
        
    for seed = 1:nr_seeds

        mvpa_dat = fmri_data_objs{seed};
        mvpa_dats{seed} = mvpa_dat;

            fprintf('\n\n');
            printhdr(['SEED #', num2str(seed), ': ', upper(seednames{seed})]);
            fprintf('\n\n');

        % DATA VISUALIZATION
        % ------------------

        fprintf('\n\n');
        printhdr('Plotting X (brain) and Y (behavioural outcome) data');
        fprintf('\n\n');

            % CON IMAGES

            h1=figure;
            
                q = zeros(size(mvpa_dat.dat,2),3);

                for subj = 1:size(mvpa_dat.dat,2)
                    this_subj_dat = mvpa_dat.dat(:,subj);
                    q(subj,:) = quantile(this_subj_dat(:),[0.025,0.5,0.975]);
                    mu = mean(mean(this_subj_dat(:)));
                    sd = std(this_subj_dat(:));
                    h1 = plot([mu-sd, mu+sd],[subj,subj],'-');
                    hold on;
                    h2 = plot(mu,subj,'o');
                    h2.Color = h1.Color;
                end

            box off
            title(['Distribution of con weights for ' fmri_data_objs{seed}.Y_names]);
            xlabel('\beta');
            ylabel('Subject');
            hold off

            p = get(gcf,'Position');
            set(gcf,'Position',[p(1:2),1024,2048],'WindowState','Maximized');
            drawnow, snapnow;

            clear subj

            % BEHAVIORAL OUTCOME

            if seed == 1

                b1=figure;

                hold off;
                b1=histogram(mvpa_dat.Y);
                box off
                title(['Histogram of ' fmri_data_objs{seed}.Y_names]);
                xlabel(fmri_data_objs{seed}.Y_names);
                ylabel('n(observations)');
                set(gcf,'WindowState','Maximized');
                drawnow, snapnow;

            end

        % RUN MODEL
        % ---------
            
            switch ml_method
                
                case 'predict'
                    
                    % CROSS-VALIDATION FOLD SELECTION

                    fprintf('\n\n');
                    printhdr('Cross-validation fold selection');
                    fprintf('\n\n');

                        switch holdout_set_method

                            case 'no_group'

                                cv = cvpartition(size(mvpa_dat.dat,2),'KFold',nfolds);
                                fold_labels = zeros(size(mvpa_dat.dat,2),1);
                                    for subj = 1:cv.NumTestSets
                                        fold_labels(cv.test(subj)) = subj;
                                    end
                                clear subj

                            case 'group'

                                group = behav_dat.(behav_group_varname);

                                cv = cvpartition(group, 'KFold',nfolds);
                                    fold_labels = zeros(size(mvpa_dat.dat,2),1);
                                    for subj = 1:cv.NumTestSets
                                        fold_labels(cv.test(subj)) = subj;
                                    end
                                clear subj

                        end % switch holdout set method
                        
                    % FIT MODEL

                    fprintf('\n\n');
                    printhdr('Fit MVPA regression model');
                    fprintf('\n\n');

                    t0 = tic;

                    [mvpa_cverr, mvpa_stats, mvpa_optout] = predict(mvpa_dat, 'algorithm_name', algorithm_name, ...
                                'nfolds', fold_labels, 'error_type', 'mse', 'parallel', 'verbose', 0);

                    t_end = toc(t0); 

                    mvpa_stats.Y_names = mvpa_dat.Y_names;
                    
                    % STORE RESULTS IN CELL ARRAY
                    
                    mvpa_stats_results{seed} = mvpa_stats;
                    
                    
                case 'oofmridataobj'
                    
                    t0 = tic;
                    
                    % DEFINE ALGORITHM

                        switch algorithm_name

                            case 'cv_pls'
                                alg = plsRegressor(); % intiate alg as a plsRegressor estimator object, other estimators in Github repo/estimators; ; numcomponents is arbitrary, but typically low for pls - we will optimize this hyperparm later

                            case 'cv_pcr'
                                alg = pcrRegressor();

                            otherwise
                                error('\nchoice of algorithm "%s" in algorithm_name option variable is not compatible with choice of machine learning method "%s" in ml_method option variable,\n Either chance method to "predict" or change algorithm to "cv_pcr" or "cv_pls"\n\n', algorithm_name, ml_method);

                        end

                        alg.fit(mvpa_dat.dat', mvpa_dat.Y); % fit alg with brain data as predictor, Y as outcome; note that fields of alg get filled
                        % NOTE: fit is not strictly necessary at this stage, but a good test
                        alg_params = alg.get_params; % get to know the hyperparams for this algorithm, which we want to optimize

                        % DEFINE FEATURE EXTRACTOR

                        featConstructor_han = @(X)([]);
                        extractVxl = fmri2VxlFeatTransformer('metadataConstructor_funhan',featConstructor_han); % initiate extractVxl as an empty fmri2VxlFeatTransformer object; other transformers in Github repo/transformers
                        extractVxl.fit(mvpa_dat); % transformer takes mvpa_data_st object as input and stores its metadata in the brainmodel property (in the .volInfo field, nifti header style data)
                        % NOTE: fit is not strictly necessary at this stage, but a good test

                        % DEFINE PIPELINE

                        fmri_pipeline = pipeline({{'featExt',extractVxl},{'alg',alg}}); % define fmri_pcr as a pipeline object including the feature transformer and the algorithm defined above; names are arbitrary
                        fmri_pipeline.fit(mvpa_dat,mvpa_dat.Y);
                        % NOTE: fit is not strictly necessary at this stage, but a good test

                        % INNER CROSS-VALIDATION FUNCTION

                        switch holdout_set_method

                            case 'group'
                                innercv = @(X,Y) cvpartition2(X.metadata_table.(behav_group_varname), 'GroupKFold', nfolds); %, 'Group', X.metadata_table.(subj_identifier_dat_st));

                            case 'onesample'
                                innercv = @(X,Y) cvpartition2(size(Y,1), 'GroupKFold', nfolds); % 'Group', X.metadata_table.(subj_identifier_dat_st)); % define innercv as handle for anonymous function cvpartition2; other partitioners in Github repo/partitioners

                        end
                        % NOTE: we use metadata_table here, since the input to bo.fit is the
                        % fmri_data_st object mvpa_dat, which has a metadata_table field

                        cv_folds = innercv(mvpa_dat,mvpa_dat.Y); 
                        % NOTE: get cross-validation folds, not strictly necessary, but a good test

                        % DEFINE BAYESIAN OPTIMIZATION

                        switch algorithm_name

                            case 'cv_pls'
                                dims = optimizableVariable('alg__numcomponents',[1,30],'Type','integer','Transform','log');

                            case 'cv_pcr'
                                dims = optimizableVariable('alg__numcomponents',[1,floor(rank(mvpa_dat.dat)*(nfolds-1)/nfolds)],'Type','integer');

                        end

                        % NOTE: Type and number of hyperparams to optimize depends on algorithm (check alg.get_params above), as well as other settings

                        bayesOptParams = {dims, 'AcquisitionFunctionName','expected-improvement-plus',...
                            'MaxObjectiveEvaluations',30, 'UseParallel', false, 'verbose',1, 'PlotFcn', {}};

                        bo = bayesOptCV(fmri_pipeline,innercv,@get_mse,bayesOptParams);
                        bo.fit(mvpa_dat,mvpa_dat.Y);
                        bo_numcomponents = bo.estimator.estimator.numcomponents;

                        % OUTER CROSS-VALIDATION FUNCTION

                        switch holdout_set_method

                            case 'group'
                                outercv = @(X,Y) cvpartition2(X.metadata_table.(behav_group_varname), 'GroupKFold', nfolds); %, 'Group', X.metadata_table.(subj_identifier_dat_st));

                            case 'onesample'
                                outercv = @(X,Y) cvpartition2(size(Y,1), 'GroupKFold', nfolds); %, 'Group', X.metadata_table.(subj_identifier_dat_st)); % define innercv as handle for anonymous function cvpartition2; other partitioners in Github repo/partitioners

                        end
                        % NOTE: we use metadata_table here, since the input to cvGS.do is the
                        % fmri_data_st object mvpa_dat, which has a metadata_table field

                        % ESTIMATE CROSS-VALIDATED MODEL PERFORMANCE

                        cvGS = crossValScore(bo, outercv, @get_mse, 'n_parallel', nfolds, 'verbose', true);
                        % NOTE: Bogdan advises not parallizing too much for the purpose of Bayesian
                        % model optimization, since each step learns from the previous one, so we
                        % only parallelize the outer cv loop with 1 core per outer cv fold,
                        % resulting in very acceptable runtimes (on LaBGAS server with 128 GB RAM)

                        cvGS.do(mvpa_dat, mvpa_dat.Y);
                        cvGS.do_null(); % fits null model - intercept only
                        fold_labels = cvGS.fold_lbls;

                        % CREATE AN FMRI_DATA OBJECT WITH THE BETAS FOR VISUALIZATION PURPOSES

                        weight_obj = bo.estimator.transformers{1}.brainModel; % empty .dat at this stage
                        weight_obj.dat = bo.estimator.estimator.B(:); % fills mdl.dat with betas
                        
                        % STORE RESULTS IN CELL ARRAY
                    
                        mvpa_stats_results{seed} = cvGS;
                        mvpa_stats_weight_objs{seed} = weight_obj;
                        
                        t_end = toc(t0); 

                    otherwise

                        error('\ninvalid option "%s" defined in ml_method_mvpa_reg_st variable, choose between "oofmridataobj" and "predict"\n\n',ml_method_mvpa_reg_st);
                    
            end


        % VISUALIZE UNTHRESHOLDED RESULTS
        % -------------------------------

        fprintf('\n\n');
        printhdr('Plotting MVPA regression results');
        fprintf('\n\n');

            % PLOT OBSERVED VERSUS PREDICTED

            fprintf('\nPLOTTING OBSERVED VERSUS PREDICTED\n');
            
            switch ml_method
                
                case 'predict'

                    fprintf('\n%s r = %0.3f\n\n', algorithm_name, corr(mvpa_stats.yfit, mvpa_dat.Y));

                    observed = mvpa_dat.Y;
                    predicted = mvpa_stats.yfit;
                    tbl = table(observed, predicted);
                    mdl = fitlm(tbl,'predicted ~ observed','RobustOpts','on');

                    figure

                    plot(mdl);
                    xlabel({['Observed ' mvpa_stats.Y_names]}); ylabel({['Estimated ' mvpa_stats.Y_names],'(cross validated)'})

                    set(gcf,'WindowState','Maximized');
                    drawnow, snapnow;
                    
                case 'oofmridataobj'
                    
                    r = zeros(size(cvGS.Y,2),1);
                    pval = zeros(size(cvGS.Y,2),1);
                    
                    for y = 1:size(cvGS.Y,2)
                        [r(y,1),pval(y,1)] = corr(cvGS.yfit{y},cvGS.Y{y});
                    end

                    cv_r = mean(r); % average r(yfit,y) over folds = cv correlation
                    clear p
                    for p = 1:size(pval,1) % set pvalue for negative correlations to 1
                        if r(p,1) < 0
                            pval(p,1) = 1;
                        end
                    end
                    cv_pval = mean(pval);

                    cv_mse = mean(cvGS.scores); % average MSE over folds = cv model performance

                    fprintf('\n%s cross-validated r = %0.3f\n\n', algorithm_name, cv_r);
                    fprintf('\n%s cross-validated p = %0.3f\n\n', algorithm_name, cv_pval);
                    fprintf('\n%s cross-validated mse = %0.3f\n\n', algorithm_name, cv_mse);

                    f1 = cvGS.plot; % plots predicted versus observed

                    set(gcf,'WindowState','Maximized');
                    drawnow, snapnow;
                    
            end

            % PLOT MONTAGE OF UNTHRESHOLDED WEIGHTS

            fprintf('\nPLOTTTING UNTHRESHOLDED WEIGHT MAPS\n');

            whmontage = 5;

            fprintf ('\nSHOWING UNTHRESHOLDED %s RESULTS, EFFECT: %s, MASK: %s\n\n', upper(algorithm_name), behav_outcome_varname, mask_string);

            figure

            o2 = canlab_results_fmridisplay([], 'compact', 'overlay', 'mni_icbm152_t1_tal_nlin_sym_09a_brainonly.img');
            
            switch ml_method
                
                case 'predict'

                    w = mvpa_stats.weight_obj;
                    
                case 'oofmridataobj'
                    
                    w = weight_obj;
                    
            end

            wreg = region(w);

            o2 = addblobs(o2, wreg, 'splitcolor',{[.1 .8 .8] [.1 .1 .8] [.9 .4 0] [1 1 0]});
            o2 = legend(o2);
            o2 = title_montage(o2, whmontage, [algorithm_name ' unthresholded ' behav_outcome_varname ' ' mask_string]);

            figtitle = sprintf('%s_unthresholded_montage_%s', algorithm_name, mask_string);
            set(gcf, 'Tag', figtitle, 'WindowState','maximized');
            drawnow, snapnow;

            clear w, clear o2, clear figtitle

    end
    
    
%% SAVE RESULTS
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('SAVING MVPA RESULTS');
fprintf('\n\n');

savefilenamedata = fullfile(resultsdir, ['mvpa_conn_stats_and_maps_', algorithm_name, '_', behav_outcome_varname, '_', maskname_short, '.mat']);
save(savefilenamedata, 'mvpa_stats_results', 'mvpa_stats_weight_objs', 'mvpa_dats','-v7.3');
fprintf('\nSaved mvpa_conn_stats_and_maps\n');

fprintf('\nFilename: %s\n', savefilenamedata);
    
    