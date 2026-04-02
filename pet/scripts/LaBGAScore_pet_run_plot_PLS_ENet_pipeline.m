%% LaBGAScore_pet_run_plot_PLS_ENet_pipeline.m
%
%
% USAGE
%
% This script serves as a simple wrapper to run the PLS-DA and Elastic Net
% neuroimaging pipeline functions and their plotting functions on parcel-wise 
% (whole-brain or ROI) PET data.
%
% The script takes excel outputs from the LaBGAScore_pet_model_TSPO_DPA714.m script
% as input, but residualizing for genotype (and if needed other covariates)
% needs to be done PRIOR to using this script!
%
% For more info, type the following in your Matlab command window
% 
% help PLSDA_neuroimaging_pipeline
% help plot_PLSDA_diagnostics_neuroimaging
% help ENet_neuroimaging_pipeline
% help plot_ENet_diagnostics_neuroimaging
% 
% or check the READMEs in the LaBGAScore Github repo
%
% https://github.com/labgas/LaBGAScore/blob/main/secondlevel/README_PLSDA_neuroimaging_pipeline.md
% https://github.com/labgas/LaBGAScore/blob/main/secondlevel/README_PLSDA_plotting.md
% https://github.com/labgas/LaBGAScore/blob/main/secondlevel/README_ENet_neuroimaging_pipeline.md
% https://github.com/labgas/LaBGAScore/blob/main/secondlevel/README_ENet_plotting.md
%
%__________________________________________________________________________
%
% Author: Lukas Van Oudenhove
% date: KU Leuven, March 2026    
%
%__________________________________________________________________________
% @(#)% LaBGAScore_pet_run_PLS_ENet_pipeline.m          v1.0
% last modified: 2026/03/07


%% ========================================================================
% 0. USER SETTINGS — EDIT THESE
% =========================================================================

% CHOOSE PIPELINE(S)

do_pls = true;
do_enet = false;


% INPUT DIRECTORIES

LaBGAScore_pet_s3_a_set_up_paths_always_run_first;

roi_resultsdir = fullfile(resultsdir,'results_master_files');

data2analyze = {'K1_ROI','K1_parcel','DV_ROI','DV_parcel'};

group_ID = 'patient';


% INPUT DATA

% pre-allocate

varnames = cell(1, size(data2analyze,2));
X_vars = cell(1, size(data2analyze,2));

% read results file

input_data = readtable(fullfile(roi_resultsdir,'Discoverie_PET_final_datafile.xlsx'));
K1_data = input_data.Properties.VariableNames(contains(input_data.Properties.VariableNames,'K1'));
DV_data = input_data.Properties.VariableNames(contains(input_data.Properties.VariableNames,'DV'));

% create cell arrays with var names for each set of input vars

varnames{1} = K1_data(1,end-13:end); % exclude parcels, only ROIs
varnames{2} = K1_data(1,1:end-14); % exclude ROIs, only parcels
varnames{3} = DV_data(1,end-13:end); % exclude parcels, only ROIs
varnames{4} = DV_data(1,1:end-14); % exclude ROIs, only parcels

% create cell arrays with X vars, and define single Y var

for x = 1:size(X_vars,2)
    X_vars{x} = table2array(input_data(:,varnames{x}));
end

Y_var = input_data.(group_ID);


% SET OPTIONS FOR PLS AND ENET PIPELINES

% Partial Least Squares
% help PLS_neuroimaging_pipeline for details

opts_PLS.outerK = 4;
opts_PLS.innerK = 4;
opts_PLS.nrepeats = 50;
opts_PLS.maxLV = 3;
opts_PLS.nPerm = 5000;
opts_PLS.nBoot = 5000;
opts_PLS.learningSteps = 6;

% Elastic Net
% help ENet_neuroimaging_pipeline for details

opts_ENet.outerK = 4;
opts_ENet.innerK = 4;
opts_ENet.nrepeats = 50;
opts_ENet.alphaGrid = [0.05 0.1 0.25 0.5 0.75 0.9 1];
opts_ENet.lambdaGrid = logspace(-3,1,25);
opts_ENet.nPerm = 5000;
opts_ENet.nBoot = 5000;
opts_ENet.learningSteps = 6;


% LOAD ATLAS

% Parcel: (selected) parcels from canlab2024

atlas = load_atlas('canlab2024');
atlas = downsample_parcellation(atlas,'labels_2');
atlas.probability_maps = [];
idx_excluded_regions = [180 181 182 183 190 191 193 201 203 209 211 214 215 216 217 218 219 220 221 222 223 229 232 245 246]; % hardcoded
atlas_reduced = atlas.remove_atlas_region(idx_excluded_regions);
atlas_reduced.fullpath = fullfile(maskdir,'canlab2024_221parcels.nii');
atlas_reduced.write('overwrite');
parcelatlasFile = fullfile(maskdir,'canlab2024_221parcels.nii');
parcelNames = atlas_reduced.labels';

% ROI: atlas object made with LaBGAScore_atlas_rois_from_atlas script

load(fullfile(maskdir,"pet_trc-DPA714_combinedTSPO.mat"));
roiatlasFile = fullfile(maskdir,'combined_TSPO.nii');
roiNames = roi_atlas.labels';


% T1 UNDERLAY FOR PLOTTING

T1 = which('fmriprep20_template.nii');
T1_obj = fmri_data(T1);
T1_obj_resample = resample_space(T1_obj,roi_atlas);
T1_obj_resample.write('fname',fullfile(maskdir,'fmriprep20_template_downsample.nii'),'overwrite');
T1_downsample = fullfile(maskdir,'fmriprep20_template_downsample.nii');


% OUTPUT DIRECTORY

pipeline_resultsdir = fullfile(resultsdir,'pls_enet_pipeline');

    if ~exist(pipeline_resultsdir,'dir')
        mkdir(pipeline_resultsdir);
    end


%% ========================================================================
% 1. CALL PIPELINE AND PLOTTING FUNCTIONS
% =========================================================================

% PLS

if do_pls
    
    PLS_results = cell(1, size(data2analyze,2));
    PLS_tables = cell(1, size(data2analyze,2));
    
    for d = 1:size(data2analyze,2)
        
        pipeline_resultssubdir = fullfile(pipeline_resultsdir,data2analyze{d});
        
            if ~exist(pipeline_resultssubdir,'dir')
                mkdir(pipeline_resultssubdir);
            end
            
        cd(pipeline_resultssubdir);

        PLS_results{d} = PLSDA_neuroimaging_pipeline(X_vars{d},Y_var,opts_PLS);
        
        [max_varY, idx_LV] = max(PLS_results{d}.varExplainedY); 
        
        fprintf('\nPlotting latent variable %d explaining %.2f%% of the variance in Y\n\n', idx_LV, max_varY*100);
        
        if contains(data2analyze{d},'ROI')
    
            PLS_tables{d} = plot_PLSDA_diagnostics_neuroimaging(PLS_results{d}, [], roiNames, roiatlasFile, ...
                'LV',idx_LV,'TopN',min(size(X_vars{d},2),20),'TopK',min(size(X_vars{d},2),20),'VIP_thresh',0.8,'stab_thresh',1.5,'MapPrctile',70,'OutPrefix',[data2analyze{d} '_PLS'],'RelaxIfEmpty',false,'UnderlayFile',T1_downsample);
            
            save_all_open_figures_smart(pipeline_resultssubdir,[data2analyze{d} '_PLS'],{'fig','svg'},true);
        
        else
            
            PLS_tables{d} = plot_PLSDA_diagnostics_neuroimaging(PLS_results{d}, [], parcelNames, parcelatlasFile, ...
                'LV',idx_LV,'TopN',min(size(X_vars{d},2),20),'TopK',min(size(X_vars{d},2),20),'VIP_thresh',0.8,'stab_thresh',1.5,'MapPrctile',70,'OutPrefix',[data2analyze{d} '_PLS'],'RelaxIfEmpty',false,'UnderlayFile',T1_downsample);
           
            save_all_open_figures_smart(pipeline_resultssubdir,[data2analyze{d} '_PLS'],{'fig','svg'},true);
            
        end
        
        clear pipeline_resultssubdir
        
    end
    
    saveplsfilename = fullfile(pipeline_resultsdir,'PLS_DA.mat');
    save(saveplsfilename, 'data2analyze','varnames','PLS_results','PLS_tables','-v7.3');

end

% ELASTIC NET

if do_enet
    
    ENet_results = cell(1, size(data2analyze,2));
    ENet_tables = cell(1, size(data2analyze,2));
    
    for d = 1:size(data2analyze,2)
        
        pipeline_resultssubdir = fullfile(pipeline_resultsdir,data2analyze{d});
        
            if ~exist(pipeline_resultssubdir,'dir')
                mkdir(pipeline_resultssubdir);
            end
            
        cd(pipeline_resultssubdir);

        ENet_results{d} = ENet_neuroimaging_pipeline(X_vars{d},Y_var,opts_ENet);
        
        if contains(data2analyze{d},'ROI')
    
            ENet_tables{d} = plot_ENet_diagnostics_neuroimaging(ENet_results{d}, X_vars{d},Y_var, roiNames, roiatlasFile, ...
                'TopN',min(size(X_vars{d},2),20),'TopK',min(size(X_vars{d},2),20),'FreqThresh',0.5,'WeightThresh',0,'MapPrctile',70,'DoPostSelection',true,'OutPrefix',[data2analyze{d} '_ENet'],'RelaxIfEmpty',false);
            
        else
            
            ENet_tables{d} = plot_ENet_diagnostics_neuroimaging(ENet_results{d}, X_vars{d},Y_var, parcelNames, parcelatlasFile, ...
                'TopN',min(size(X_vars{d},2),20),'TopK',min(size(X_vars{d},2),20),'FreqThresh',0.5,'WeightThresh',0,'MapPrctile',70,'DoPostSelection',true,'OutPrefix',[data2analyze{d} '_ENet'],'RelaxIfEmpty',false);
            
        end
        
        clear pipeline_resultssubdir
    
    end
    
    saveenetfilename = fullfile(pipeline_resultsdir,'ENet.mat');
    save(saveenetfilename, 'data2analyze','varnames','ENet_results','ENet_tables','-v7.3');

end

cd(rootdir)


%% ========================================================================
% HELPER FUNCTION
% =========================================================================

function outFiles = save_all_open_figures_smart(outDir, prefix, formats, closeFigs)
% save_all_open_figures_smart
%
% Save all open MATLAB figures with smart filenames.
%
% Priority for naming:
%   1) figure Name
%   2) axes title
%   3) fallback: fig01, fig02, ...
%
% INPUTS
%   outDir     - output directory (default: pwd)
%   prefix     - filename prefix (default: 'Figure')
%   formats    - cell array, e.g. {'fig','svg'} (default: {'fig','svg'})
%   closeFigs  - logical (default: false)
%
% OUTPUT
%   outFiles   - struct with saved file paths
%
% EXAMPLE
%   save_all_open_figures_smart('results','PLSDA',{'fig','svg'},false)

if nargin < 1 || isempty(outDir)
    outDir = pwd;
end
if nargin < 2 || isempty(prefix)
    prefix = 'Figure';
end
if nargin < 3 || isempty(formats)
    formats = {'fig','svg'};
end
if nargin < 4
    closeFigs = false;
end

if ischar(formats) || isstring(formats)
    formats = cellstr(formats);
end

if ~exist(outDir,'dir')
    mkdir(outDir);
end

figs = findall(0,'Type','figure');

if isempty(figs)
    warning('No open figures to save.');
    outFiles = struct([]);
    return
end

% Sort by figure number
[~,idx] = sort([figs.Number]);
figs = figs(idx);

outFiles = struct('figureNumber',{},'name',{},'files',{});
usedNames = containers.Map('KeyType','char','ValueType','double');

for i = 1:numel(figs)

    h = figs(i);

    % ---------- Get smart name ----------
    name = get(h,'Name');

    if isempty(name)
        name = get_axes_title(h);
    end

    if isempty(name)
        name = sprintf('fig%02d',i);
    end

    name = sanitize_filename(name);
    prefixSafe = sanitize_filename(prefix);

    % Avoid duplicates
    if isKey(usedNames,name)
        usedNames(name) = usedNames(name) + 1;
        name = sprintf('%s_%02d', name, usedNames(name));
    else
        usedNames(name) = 1;
    end

    base = fullfile(outDir, sprintf('%s_%s', prefixSafe, name));

    saved = {};

    % ---------- Save formats ----------
    for f = 1:numel(formats)

        fmt = lower(formats{f});

        switch fmt
            case 'fig'
                file = [base '.fig'];
                savefig(h,file,'compact');

            case 'svg'
                file = [base '.svg'];
                try
                    exportgraphics(h,file,'ContentType','vector');
                catch
                    print(h,file,'-dsvg');
                end

            case 'png'
                file = [base '.png'];
                exportgraphics(h,file,'Resolution',300);

            case 'pdf'
                file = [base '.pdf'];
                exportgraphics(h,file,'ContentType','vector');

            otherwise
                warning('Unsupported format: %s',fmt);
                continue
        end

        saved{end+1} = file;
    end

    outFiles(i).figureNumber = h.Number;
    outFiles(i).name = name;
    outFiles(i).files = saved;

    if closeFigs
        close(h);
    end
end

fprintf('Saved %d figure(s) to %s\n', numel(figs), outDir);

end


% =========================
% Helper: get axes title
% =========================
function t = get_axes_title(hFig)

t = '';

ax = findall(hFig,'Type','axes');

for k = 1:numel(ax)
    try
        titleStr = get(get(ax(k),'Title'),'String');
    catch
        titleStr = '';
    end

    if iscell(titleStr)
        titleStr = strjoin(titleStr,' ');
    end

    if isstring(titleStr)
        titleStr = char(titleStr);
    end

    titleStr = strtrim(titleStr);

    if ~isempty(titleStr)
        t = titleStr;
        return
    end
end
end


% =========================
% Helper: sanitize filename
% =========================
function s = sanitize_filename(s)

if isstring(s)
    s = char(s);
end

s = strtrim(s);

% Replace spaces
s = regexprep(s,'\s+','_');

% Remove problematic chars
s = regexprep(s,'[^\w\-\.]','_');

% Collapse underscores
s = regexprep(s,'_+','_');

% Trim
s = regexprep(s,'^_+|_+$','');

if isempty(s)
    s = 'figure';
end

end