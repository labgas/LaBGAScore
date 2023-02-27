%%% LaBGAScore_secondlevel_ooFmriDataObjML_example.m
%
% This is merely and example script to illustrate the object-oriented ML
% toolkit for fmri_data objects (and beyond), developed by Bogdan Petre
% @CANlab.
%
% This is a toolkit inspired by Python's scikit-learn, including more
% flexible options for algorithm and feature selection, 
% hyperparameter optimization, nested cross-validation, etc, compared to CANlab's predict function. 
% However, it does require more advanced programming skills and understanding the logic
% of the method.
%
% Dependency: https://github.com/canlab/ooFmriDataObjML
% Also lots of useful examples there!
%
% Tutorial: https://canlab.github.io/_pages/canlab_pipelines_walkthrough/estimateBestRegionPerformance.html
% Method also implemented in
% https://github.com/labgas/CANlab_help_examples/blob/master/Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/c2f_run_MVPA_regression_single_trial.m,
% including some more info as well
% 
%__________________________________________________________________________
%
% author: lukas.vanoudenhove@kuleuven.be & bogpetre@gmail.com
% date:   June, 2022
%__________________________________________________________________________
% @(#)% LaBGAScore_secondlevel_ooFmriDataObjML_example     v2.0        
% last modified: 2022/08/04


%% EXAMPLE 1: OPTIMIZATION USING GRID SEARCH
%-------------------------------------------

% DEFINE ALGORITHM
%-----------------
alg = plsRegressor('numcomponents',5); % intiate alg as a plsRegressor estimator object, other estimators in Github repo/estimators; ; numcomponents is arbitrary, but typically low for pls - we will optimize this hyperparm later
alg.fit(fmri_dat.dat', fmri_dat.Y); % fit alg with brain data as predictor, Y as outcome; note that fields of alg get filled
% NOTE: fit is not strictly necessary at this stage, but a good test
alg_params = alg.get_params; % get to know the hyperparams for this algorithm, which we want to optimize

% DEFINE FEATURE EXTRACTOR
%-------------------------
extractVxl = fmri2VxlFeatTransformer; % initiate extractVxl as an empty fmri2VxlFeatTransformer object; other transformers in Github repo/transformers
extractVxl.fit(fmri_dat); % transformer takes fmri_data_st object as input and stores its metadata in the brainmodel property (in the .volInfo field, nifti header style data)
% NOTE: fit is not strictly necessary at this stage, but a good test

% DEFINE PIPELINE
%----------------
fmri_pipeline = pipeline({{'featExt',extractVxl},{'alg',alg}}); % define fmri_pcr as a pipeline object including the feature transformer and the algorithm defined above; names are arbitrary
fmri_pipeline.fit(fmri_dat,fmri_dat.Y);
% NOTE: fit is not strictly necessary at this stage, but a good test

% INNER CROSS-VALIDATION FUNCTION
%--------------------------------
switch holdout_set_method_mvpa_reg_st

    case 'group'
        innercv = @(X,Y) cvpartition2(X.metadata_table.(group_identifier), 'GroupKFold', nfolds_mvpa_reg_st, 'Group', X.metadata_table.subject_id);

    case 'onesample'
        innercv = @(X,Y) cvpartition2(size(Y,1), 'GroupKFold', nfolds_mvpa_reg_st, 'Group', X.metadata_table.subject_id); % define innercv as handle for anonymous function cvpartition2; other partitioners in Github repo/partitioners

end
% NOTE: we use metadata_table here, since the input to go.fit is the
% fmri_data_st object fmri_dat, which has a metadata_table field

cv_folds = innercv(fmri_dat,fmri_dat.Y); 
% NOTE: get cross-validation folds, not strictly necessary, but a good test

% DEFINE OPTIMIZATION GRID
%-------------------------
gridPoints = table([2 4 5 9 20]','VariableNames',{'alg__numcomponents'}); 
% NOTES 
% 1. the double underscore linking name of algorithm in pipeline with its parameter we want to optimize
% 2. the choice of the discrete points over which we want to perform our
% grid search is based on knowledge of the algorithm, in this case PLS
% which typically has a rather low number of components, with the first
% ones explaining most of the variance

go = gridSearchCV(fmri_pipeline,gridPoints,innercv,@get_mse,'verbose',true); % see help gridSearchCV for inputs; scorers in Github repo/scorer; other optimizers in Github repo/estimators
go.fit(fmri_dat,fmri_dat.Y);

% OUTER CROSS-VALIDATION FUNCTION
%--------------------------------
switch holdout_set_method_mvpa_reg_st

    case 'group'
        outercv = @(X,Y) cvpartition2(X.metadata_table.(group_identifier), 'GroupKFold', nfolds_mvpa_reg_st, 'Group', X.metadata_table.subject_id);

    case 'onesample'
        outercv = @(X,Y) cvpartition2(size(Y,1), 'GroupKFold', nfolds_mvpa_reg_st, 'Group', X.metadata_table.subject_id); % define innercv as handle for anonymous function cvpartition2; other partitioners in Github repo/partitioners

end
% NOTE: we use metadata_table here, since the input to cvGS.do is the
% fmri_data_st object fmri_dat, which has a metadata_table field

% ESTIMATE CROSS-VALIDATED MODEL PERFORMANCE
%-------------------------------------------
cvGS = crossValScore(go, outercv, @get_mse, 'verbose', true); % see help/methods crossValScore; more crossvalidators in Github_repo/crossValidators - check out crossValPredict to get cross-validated prediction estimates at this stage
cvGS.do(fmri_dat, fmri_dat.Y);
cvGS.do_null(); % fits null model - intercept only

figure
cvGS.plot; % plots predicted versus observed
% average MSE over folds = cv model performance

% CREATE AN FMRI_DATA OBJECT WITH THE BETAS FOR VISUALIZATION PURPOSES
%---------------------------------------------------------------------
mdl = go.estimator.transformers{1}.brainModel; % empty .dat at this stage
mdl.dat = go.estimator.estimator.B(:); % fills mdl.dat with betas


%% EXAMPLE 2: BAYESIAN OPTIMIZATION & HARMONIZATION
%--------------------------------------------------------------------------
% NOTE: 
% This example is not only an illustration of how to use Bayesian rather than
% gridsearch-based optimization (which is the preferred option), but also
% of an alternative order of the different steps, which is a bit more
% complicated/convoluted and may be less clean, but could be faster. It
% also illustrates how to build in z-scoring voxels and/or images in your
% pipeline

% DEFINE ALGORITHM
%-----------------
alg = plsRegressor();
alg.fit(fmri_dat.dat', fmri_dat.Y);
% NOTE: fit is not strictly necessary at this stage, but a good test
alg_params = alg.get_params; % get to know the hyperparams for this algorithm, which we want to optimize

% INNER CROSS-VALIDATION FUNCTION
%--------------------------------
switch holdout_set_method_mvpa_reg_st

    case 'group'
        innercv = @(X,Y) cvpartition2(X.metadata.(group_identifier), 'GroupKFold', nfolds_mvpa_reg_st, 'Group', X.metadata.(subj_identifier));

    case 'onesample'
        innercv = @(X,Y) cvpartition2(size(Y,1), 'GroupKFold', nfolds_mvpa_reg_st, 'Group', X.metadata.(subj_identifier)); % define innercv as handle for anonymous function cvpartition2; other partitioners in Github repo/partitioners

end
% NOTE: we use metadata.() here rather than metadata_table.() here because bo is downstream of
% a feature extractor in our pipeline below

% DEFINE BAYESIAN OPTIMIZATION
%-----------------------------
dims = optimizableVariable('numcomponents',[1,30],'Type','integer','Transform','log');
% NOTE: Type and number of hyperparams to optimize depends on algorithm (check alg.get_params above), as well as other settings
% Example for PCR: up to 75% of width of fmri_dat.dat, and no log transform of search space

bayesOptParams = {dims, 'AcquisitionFunctionName','expected-improvement-plus',...
    'MaxObjectiveEvaluations',30, 'UseParallel', false, 'verbose',1, 'PlotFcn', {}};
bo = bayesOptCV(alg,innercv,@get_mse,bayesOptParams);

% DEFINE OTHER COMPONENTS OF PIPELINE
%------------------------------------
switch holdout_set_method_mvpa_reg_st

    case 'group'
        featConstructor_han = @(X) table(X.metadata_table.(subj_identifier),X.metadata_table.(group_identifier),'VariableNames',{subj_identifier,group_identifier});
        
    case 'onesample'
        featConstructor_han = @(X) table(X.metadata_table.(subj_identifier),'VariableNames',{subj_identifier});
        
end

extractVxl = fmri2VxlFeatTransformer('metadataConstructor_funhan',featConstructor_han);

zTransVxl = zscoreVxlTransformer(@(X) X.metadata_table.(subj_identifier));
zTransImg = functionTransformer(@(x1) rescale(x1,'zscoreimages'));

% DEFINE PIPELINE
%------------------------------------
% fmri_pipeline = pipeline({{'zscorevxl', zTransVxl},{'zscoreimg', zTransImg},...
%                     {'featExt',extractVxl},{'bo_pls',bo}});
% NOTE: this is a pipeline including zscoring voxels and images
fmri_pipeline = pipeline({{'featExt',extractVxl},{'bo_pls',bo}});
% feat = features(fmri_dat.dat', fmri_dat.metadata_table.subject_id); 
% NOTES
% 1. this is an "extended double" that is just a double with metadata in the dat.metadata field 
% 2. this is an alternative to the use of featConstructor_han in the pipeline above, which
% can go into fmri_ppipeline.fit as an argument instead of fmri_dat
fmri_pipeline.fit(fmri_dat,fmri_dat.Y);
% fmri_pipeline.fit(feat,fmri_dat.Y); % see notes above

% OUTER CROSS-VALIDATION FUNCTION
%--------------------------------
switch holdout_set_method_mvpa_reg_st

    case 'group'
        outercv = @(X,Y) cvpartition2(X.metadata_table.(group_identifier), 'GroupKFold', nfolds_mvpa_reg_st, 'Group', X.metadata_table.subject_id);

    case 'onesample'
        outercv = @(X,Y) cvpartition2(size(Y,1), 'GroupKFold', nfolds_mvpa_reg_st, 'Group', X.metadata_table.subject_id); % define innercv as handle for anonymous function cvpartition2; other partitioners in Github repo/partitioners

end
% NOTE: we use metadata_table here, since the input to cvGS.do is the
% fmri_data_st object fmri_dat, which has a metadata_table field

% ESTIMATE CROSS-VALIDATED MODEL PERFORMANCE
%-------------------------------------------
cvGS = crossValScore(fmri_pipeline, outercv, @get_mse, 'n_parallel', nfolds_mvpa_reg_st, 'verbose', true);
% NOTE: Bogdan advises not parallizing too much for the purpose of Bayesian
% model optimization, since each step learns from the previous one, so we
% only parallelize the outer cv loop with 1 core per outer cv fold,
% resulting in very acceptable runtimes (on LaBGAS server with 128 GB RAM)

cvGS.do(fmri_dat, fmri_dat.Y);
cvGS.do_null();

figure
cvGS.plot; % plots predicted versus observed
% average MSE over folds = cv model performance

% CREATE AN FMRI_DATA OBJECT WITH THE BETAS FOR VISUALIZATION PURPOSES
%---------------------------------------------------------------------
mdl = fmri_pls.transformers{end}.brainModel; % empty .dat at this stage
mdl.dat = fmri_pls.estimator.estimator.B(:); % fills mdl.dat with betas
