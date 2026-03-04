%% load data

data = readtable("Discoverie_PET_final_datafile.xlsx");
K1_vars = data.Properties.VariableNames(contains(data.Properties.VariableDescriptions,'K1'));
K1_vars = K1_vars(1,1:end-14); % exclude ROIs, only parcels
DV_vars = data.Properties.VariableNames(contains(data.Properties.VariableDescriptions,'DV'));
DV_vars = DV_vars(1,1:end-14); % exclude ROIs, only parcels

%% input for functions

K1 = table2array(data(:,K1_vars)); % X1
DV = table2array(data(:,DV_vars)); % X2
group = data.patient; % Y
n = height(data); % n

opts.outerK = 4;
opts.innerK = 4;
opts.nrepeats = 50;
opts.maxLV = 4;
opts.nPerm = 5000;
opts.nBoot = 5000;
opts.learningSteps = 6;

opts_ENet.outerK = 4;
opts_ENet.innerK = 4;
opts_ENet.nrepeats = 10;
opts_ENet.alphaGrid = [0.05 0.1 0.25 0.5 0.75 0.9 1];
opts_ENet.lambdaGrid = logspace(-4,0,20);
opts_ENet.nPerm = 100;
opts_ENet.nBoot = 100;
opts_ENet.learningSteps = 6;

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

results_K1_TSPO = TSPO_PLSDA_pipeline(K1,group,opts);
results_K1_TSPO_ENet = TSPO_ENet_pipeline(K1,group,opts_ENet);
results_K1_generic = PLSDA_neuroimaging_pipeline(K1,group,opts);
results_K1_generic_ENet = ENet_neuroimaging_pipeline(K1,group,opts_ENet);

plot_ENet_diagnostics_neuroimaging(results_K1_generic_ENet, K1, group, roiNames, atlasFile, ...
    'TopN',20,'TopK',20,'FreqThresh',0.5,'WeightThresh',0,'MapPrctile',70,'DoPostSelection',true,'OutPrefix','generic_ENet','RelaxIfEmpty',true);

ROI_table = plot_PLSDA_diagnostics_neuroimaging(results_K1_TSPO, [], roiNames, atlasFile, ...
    'LV',2,'TopN',20,'TopK',20,'VIP_thresh',0.8,'stab_thresh',1.5,'MapPrctile',70,'OutPrefix','generic_PLSDA','RelaxIfEmpty',true);
