%% EXAMPLE 1: OPTIMIZATION USING GRID SEARCH

% DEFINE ALGORITHM AND FEATURE EXTRACTOR
alg = plsRegressor('numcomponents',5); % intiate alg as a plsRegressor estimator object, other estimators in Github repo/estimators; numcomponents is arbitrary, but typically low for pls - we will optimize this hyperparm later
alg.fit(fmri_dat.dat', fmri_dat.Y); % fit alg with brain data as predictor, Y as outcome; note that fields of alg get filled
extractVxl = fmri2VxlFeatTransformer; % initiate extractVxl as an (empty?) fmri2VxlFeatTransformer object; other transformers in Github repo/transformers
extractVxl.fit(fmri_dat); % transformer takes fmri_data_st object as input and stores its metadata in the brainmodel property (in the .volInfo field, nifti header style data)

% DEFINE PIPELINE
fmri_pls = pipeline({{'featExt',extractVxl},{'pls',alg}}); % define fmri_pls as a pipeline object including the feature transformer and the algorithm defined above; names are arbitrary
fmri_pls.fit(fmri_dat,fmri_dat.Y); % fit pipeline

% INNER CROSS-VALIDATION LOOP
innercv = @(X,Y) cvpartition2(length(Y), 'GroupKFold', 2, 'Group', X.metadata_table.participant_id); % define innercv as handle for anonymous function cvpartition2; other partitioners in Github repo/partitioners
innercv(fmri_dat,fmri_dat.Y); % get cross-validation folds

% DEFINE OPTIMIZATION GRID
gridPoints = table([2 4 5 9 20]','VariableNames',{'pls__numcomponents'}); % note the double underscore linking name of algorithm in pipeline with its parameter we want to optimize
go = gridSearchCV(fmri_pls,gridPoints,innercv,@get_mse,'verbose',true); % see help gridSearchCV for inputs; scorers in Github repo/scorer; other optimizers in Github repo/estimators

go.fit(fmri_dat,fmri_dat.Y);
mdl = go.estimator.transformers{1}.brainModel; % empty .dat at this stage
mdl.dat = go.estimator.estimator.B; % fills mdl.dat with betas

% OUTER CROSS-VALIDATION LOOP
outercv = @(X,Y) cvpartition2(length(Y), 'GroupKFold', 3, 'Group', X.metadata_table.participant_id); % same as above

% ESTIMATE CROSS-VALIDATED MODEL PERFORMANCE
cvGS = crossValScore(bo, outercv, @get_mse, 'verbose', true); % see help/methods crossValScore; more crossvalidators in Github_repo/crossValidators - check out crossValPredict to get cross-validated prediction estimates at this stage
cvGS.do(fmri_dat, fmri_dat.Y);
cvGS.do_null(); % fits null model - intercept only
f1 = cvGS.plot; % plots predicted versus observed
% average MSE over folds = cv model performance


%% EXAMPLE 2: BAYESIAN OPTIMIZATION & HARMONIZATION
%--------------------------------------------------------------------------

alg = plsRegressor();
alg.fit(fmri_dat.dat', fmri_dat.Y);

innercv = @(X,Y) cvpartition2(length(Y), 'GroupKFold', 5, 'Group', X.metadata.participant_id);

dims = optimizableVariable('numcomponents',[1,30],'Type','integer','Transform','log');
bayesOptParams = {dims, 'AcquisitionFunctionName','expected-improvement-plus',...
    'MaxObjectiveEvaluations',30,'verbose',1, 'PlotFcn', {}};
bo = bayesOptCV(alg,innercv,@get_mse,bayesOptParams);


featConstructor_han = @(X) table(X.metadata_table.participant_id,'VariableNames',{'participant_id'});

extractVxl = fmri2VxlFeatTransformer('metadataConstructor_funhan',featConstructor_han);
zTransVxl = zscoreVxlTransformer(@(X) X.metadata_table.participant_id);

fmri_pls = pipeline({{'zscorevxl', zTransVxl},{'zscoreimg', functionTransformer(@(x1)rescale(x1,'zscoreimages'))},...
    {'featExt',extractVxl},{'bo_pls',bo}});
fmri_pls.fit(fmri_dat,fmri_dat.Y);

outercv = @(X,Y) cvpartition2(length(Y), 'GroupKFold', 5, 'Group', X.metadata_table.participant_id);

cvGS = crossValScore(fmri_pls, outercv, @get_mse, 'verbose', true);
cvGS.do(fmri_dat, fmri_dat.Y);
cvGS.do_null();
cvGS.plot();