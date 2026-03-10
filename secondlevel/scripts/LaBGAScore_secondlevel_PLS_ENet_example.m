%% load data

data = readtable("Discoverie_PET_final_datafile.xlsx");
K1_data = data.Properties.VariableNames(contains(data.Properties.VariableDescriptions,'K1'));
K1_parcels = K1_data(1,1:end-14); % exclude ROIs, only parcels
K1_ROIs = K1_data(1,end-13:end); % exclude parcels, only ROIs
DV_data = data.Properties.VariableNames(contains(data.Properties.VariableDescriptions,'DV'));
DV_parcels = DV_data(1,1:end-14); % exclude ROIs, only parcels
DV_ROIs = DV_data(1,end-13:end); % exclude parcels, only ROIs

data_rif = readtable("roi_means_table_rifaximin.xlsx");

data_graphvars = readtable("graphvars_GSR.xlsx");
data_graphvars_D30 = data_graphvars.Properties.VariableNames(contains(data_graphvars.Properties.VariableDescriptions,'D30'));
data_graphvars_char_path = data_graphvars.Properties.VariableNames(contains(data_graphvars.Properties.VariableDescriptions,'char_path'));
data_graphvars_eff_global = data_graphvars.Properties.VariableNames(contains(data_graphvars.Properties.VariableDescriptions,'efficiency_global'));
data_graphvars_bc = data_graphvars.Properties.VariableNames(contains(data_graphvars.Properties.VariableDescriptions,'betweenness_centrality'));
data_graphvars_D30_full = [data_graphvars_D30, data_graphvars_char_path, data_graphvars_eff_global, data_graphvars_bc];
data_graphvars_D35 = data_graphvars.Properties.VariableNames(contains(data_graphvars.Properties.VariableDescriptions,'D35'));
data_graphvars_D35_full = [data_graphvars_D35, data_graphvars_char_path, data_graphvars_eff_global, data_graphvars_bc];
data_graphvars_D40 = data_graphvars.Properties.VariableNames(contains(data_graphvars.Properties.VariableDescriptions,'D40'));
data_graphvars_D40_full = [data_graphvars_D40, data_graphvars_char_path, data_graphvars_eff_global, data_graphvars_bc];
data_graphvars_D45 = data_graphvars.Properties.VariableNames(contains(data_graphvars.Properties.VariableDescriptions,'D45'));
data_graphvars_D45_full = [data_graphvars_D45, data_graphvars_char_path, data_graphvars_eff_global, data_graphvars_bc];
data_graphvars_D50 = data_graphvars.Properties.VariableNames(contains(data_graphvars.Properties.VariableDescriptions,'D50'));
data_graphvars_D50_full = [data_graphvars_D50, data_graphvars_char_path, data_graphvars_eff_global, data_graphvars_bc];

data_graphvars_all = data_graphvars(:,4:543);

%% input for functions

K1_parcels_X = table2array(data(:,K1_parcels)); % X1
DV_parcels_X = table2array(data(:,DV_parcels)); % X2
K1_ROIs_X = table2array(data(:,K1_ROIs)); % X1
DV_ROIs_X = table2array(data(:,DV_ROIs)); % X2
group = data.patient; % Y
n = height(data); % n

ROIs_rif = table2array(data_rif(:,1:end-1));
group_rif = data_rif.group; % Y
n = height(data_rif); % n

graphvars_D30 = table2array(data_graphvars(:,data_graphvars_D30_full));
graphvars_D35 = table2array(data_graphvars(:,data_graphvars_D35_full));
graphvars_D40 = table2array(data_graphvars(:,data_graphvars_D40_full));
graphvars_D45 = table2array(data_graphvars(:,data_graphvars_D45_full));
graphvars_D50 = table2array(data_graphvars(:,data_graphvars_D50_full));
graphvars_all = table2array(data_graphvars_all);
group_graphvars = data_graphvars.Group;

opts.outerK = 4;
opts.innerK = 4;
opts.nrepeats = 50;
opts.maxLV = 4;
opts.nPerm = 10000;
opts.nBoot = 10000;
opts.learningSteps = 6;

opts_ENet.outerK = 4;
opts_ENet.innerK = 4;
opts_ENet.nrepeats = 50;
opts_ENet.alphaGrid = [0.05 0.1 0.25 0.5 0.75 0.9 1];
opts_ENet.lambdaGrid = logspace(-3,1,25);
opts_ENet.nPerm = 10000;
opts_ENet.nBoot = 10000;
opts_ENet.learningSteps = 6;

opts_rif.outerK = 8;
opts_rif.innerK = 4;
opts_rif.nrepeats = 50;
opts_rif.maxLV = 4;
opts_rif.nPerm = 10000;
opts_rif.nBoot = 10000;
opts_rif.learningSteps = 6;

opts_rif_ENet.outerK = 8;
opts_rif_ENet.innerK = 4;
opts_rif_ENet.nrepeats = 50;
opts_rif_ENet.alphaGrid = [0.05 0.1 0.25 0.5 0.75 0.9 1];
opts_rif_ENet.lambdaGrid = logspace(-3,1,25);
opts_rif_ENet.nPerm = 5000;
opts_rif_ENet.nBoot = 5000;
opts_rif_ENet.learningSteps = 6;

atlas = load_atlas('canlab2024');
atlas = downsample_parcellation(atlas,'labels_2');
atlas.probability_maps = [];
idx_excluded_regions = [180 181 182 183 190 191 193 201 203 209 211 214 215 216 217 218 219 220 221 222 223 229 232 245 246];
atlas_reduced = atlas.remove_atlas_region(idx_excluded_regions);
atlas_reduced.fullpath = 'canlab2024_221parcels.nii';
% atlas_reduced.write
atlasFile = atlas_reduced.fullpath;
roiNames = atlas_reduced.labels';

%% call function

results_K1_TSPO = TSPO_PLSDA_pipeline(K1_parcels_X,group,opts);
results_DV_TSPO = TSPO_PLSDA_pipeline(DV_parcels_X,group,opts);
results_K1_TSPO_ENet = TSPO_ENet_pipeline(K1_parcels_X,group,opts_ENet);
results_K1_ROIs = PLSDA_neuroimaging_pipeline(K1_ROIs_X,group,opts);
results_K1_ROIs_ENet = ENet_neuroimaging_pipeline(K1_ROIs_X,group,opts_ENet);
results_DV_ROIs = PLSDA_neuroimaging_pipeline(DV_ROIs_X,group,opts);
results_DV_ROIs_ENet = ENet_neuroimaging_pipeline(DV_ROIs_X,group,opts_ENet);

plot_ENet_diagnostics_neuroimaging(results_K1_TSPO_ENet, K1_parcels_X, group, roiNames, atlasFile, ...
    'TopN',20,'TopK',20,'FreqThresh',0.5,'WeightThresh',0,'MapPrctile',70,'DoPostSelection',true,'OutPrefix','generic_ENet','RelaxIfEmpty',true);

parcel_table_K1 = plot_PLSDA_diagnostics_neuroimaging(results_K1_TSPO, [], roiNames, atlasFile, ...
    'LV',2,'TopN',20,'TopK',20,'VIP_thresh',0.8,'stab_thresh',1.5,'MapPrctile',70,'OutPrefix','generic_PLSDA','RelaxIfEmpty',true);

parcel_table_DV = plot_PLSDA_diagnostics_neuroimaging(results_DV_TSPO, [], roiNames, atlasFile, ...
    'LV',1,'TopN',20,'TopK',20,'VIP_thresh',0.8,'stab_thresh',1.5,'MapPrctile',70,'OutPrefix','generic_PLSDA','RelaxIfEmpty',true);

results_rif = PLSDA_neuroimaging_pipeline(ROIs_rif,group_rif,opts_rif);
results_rif_ENet = ENet_neuroimaging_pipeline(ROIs_rif,group_rif,opts_rif_ENet);

results_graphvars = PLSDA_neuroimaging_pipeline(graphvars_D30,group_graphvars,opts_rif);
results_graphvars_ENet = ENet_neuroimaging_pipeline(graphvars_D30,group_graphvars,opts_rif_ENet);
results_graphvars_D35 = PLSDA_neuroimaging_pipeline(graphvars_D35,group_graphvars,opts_rif);
results_graphvars_D35_ENet = ENet_neuroimaging_pipeline(graphvars_D35,group_graphvars,opts_rif_ENet);
results_graphvars_D40 = PLSDA_neuroimaging_pipeline(graphvars_D40,group_graphvars,opts_rif);
results_graphvars_D40_ENet = ENet_neuroimaging_pipeline(graphvars_D40,group_graphvars,opts_rif_ENet);
results_graphvars_D45 = PLSDA_neuroimaging_pipeline(graphvars_D45,group_graphvars,opts_rif);
results_graphvars_D50 = PLSDA_neuroimaging_pipeline(graphvars_D50,group_graphvars,opts_rif);
results_graphvars_all = PLSDA_neuroimaging_pipeline(graphvars_all,group_graphvars,opts_rif);