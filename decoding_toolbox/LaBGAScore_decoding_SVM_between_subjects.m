%% LaBGAScore_decoding_SVM_between_subjects
%
% BETWEEN-SUBJECT SVM DECODING WITH PERMUTATION-BASED TFCE INFERENCE
%
% This script implements a complete second-level between-subject MVPA
% decoding pipeline using The Decoding Toolbox (TDT). It performs group-level
% SVM classification (patients vs controls), permutation-based statistical
% inference, and Threshold-Free Cluster Enhancement (TFCE) for spatially
% informed whole-brain, ROI, or searchlight decoding analyses.
%
% Statistical inference is non-parametric throughout and based on empirical
% permutation null distributions. TFCE constitutes the primary inferential
% framework, supplemented by voxel-wise permutation p-values and FDR
% correction for interpretability.
%
% The pipeline is designed for reproducibility, scalability, and integration
% with the Canlab neuroimaging framework.
%
% -------------------------------------------------------------------------
% OVERVIEW
% -------------------------------------------------------------------------
% The script executes the following steps:
%
% 1. Data Preparation
%    - Loads phenotype information (participant IDs, group labels).
%    - Constructs balanced K-fold cross-validation (CV) partitions
%      separately for patients and controls.
%    - Assigns chunk labels and class labels to each subject.
%    - Loads first-level contrast images for all subjects.
%
% 2. Decoding Configuration (TDT)
%    Between-subject linear SVM decoding using The Decoding Toolbox,
%    supporting three analysis modes:
%
%      * 'searchlight' : spherical multivoxel decoding maps
%      * 'roi'         : decoding within multiple predefined ROIs
%      * 'wholebrain'  : decoding using a full-brain mask
%
%    Key configuration choices include:
%      - Linear SVM classifier
%      - Balanced cross-validation folds
%      - Fold-wise z-scoring for real decoding (training-only scaling)
%      - Explicit handling of unbalanced group sizes
%      - Performance metric: AUC minus chance
%
% 3. Real (Unpermuted) Decoding
%    - Runs decoding once using the true group labels.
%    - Extracts the voxel-wise decoding statistic vector (real_vec),
%      consistent across analysis modes.
%
% 4. Permutation Design Generation
%    - Generates N permutation designs using TDT's permutation framework.
%    - Permutations preserve the cross-validation structure while
%      destroying the label–group relationship.
%
% 5. Fast Permutation Decoding (In-Memory)
%    - Executes permutation decodings using parfor for efficiency.
%    - Uses linear SVM decoding with fixed (global) z-scaling.
%    - Avoids writing intermediate output folders or NIfTI files.
%    - Collects the empirical null distribution in memory as:
%
%          all_perm_results [nVoxels × nPermutations]
%
%    - Live progress and ETA reporting are provided during permutation runs.
%
% 6. Voxel-wise Permutation Statistics
%    - Computes empirical voxel-wise p-values from the permutation null:
%
%          p_unc(v) = ( #{perm >= real} + 1 ) / (nPerm + 1)
%
%    - Generates uncorrected and FDR-corrected voxel-wise p-maps for
%      complementary assessment of decoding effects.
%
% 7. TFCE-Based Inference (Primary)
%    - Converts decoding statistics to voxel-wise t-maps using the
%      permutation-derived mean and standard deviation.
%    - Applies Threshold-Free Cluster Enhancement (TFCE) using pTFCE.
%    - Computes a permutation-based TFCE null distribution (in parallel).
%
%    TFCE inference includes:
%      - Global FWE-corrected TFCE p-value (max-statistic)
%      - Voxel-wise TFCE p-values (uncorrected)
%      - Voxel-wise FDR-corrected TFCE p-values
%      - Voxel-wise FWE-corrected TFCE p-values (using the TFCE max-statistic)
%
%    TFCE results constitute the primary inferential outcome of this
%    pipeline, providing spatially informed, permutation-based control
%    of family-wise error.
%
% 8. Canlab Output Generation
%    - Converts decoding statistics and p-maps into Canlab-compatible
%      representations:
%        * fmri_data_st objects
%        * region objects with atlas-based labels
%        * summary tables for reporting and visualization
%
% -------------------------------------------------------------------------
% INPUT REQUIREMENTS
% -------------------------------------------------------------------------
% - MATLAB R2018a or newer
% - The Decoding Toolbox (TDT)
% - pTFCE toolbox on MATLAB path
% - SPM12 (recommended for NIfTI I/O)
% - NIfTI utilities (load_nii / save_nii)
% - Canlab toolbox (fmri_data_st, region, atlas utilities)
%
% - phenotype.csv in BIDSdir containing:
%      * participant identifiers
%      * group labels (patient / control)
%
% - First-level contrast images (con_XXXX.nii) for all subjects
% - Valid brain mask (whole-brain / searchlight) or ROI masks
%
% -------------------------------------------------------------------------
% OUTPUTS
% -------------------------------------------------------------------------
% All outputs are written to tdt_resultsdir/, including:
%
% - NIfTI files:
%      * decoding statistic maps
%      * voxel-wise permutation p-maps
%      * TFCE-enhanced maps
%      * TFCE-based p-maps (voxel-wise, FDR, and FWE)
%
% - MATLAB files:
%      * combined_permutations.mat
%          - real_vec
%          - all_perm_results
%          - analysis metadata
%
% - Canlab objects:
%      * fmri_data_st objects (statistics and p-maps)
%      * region objects (TFCE- or voxel-wise thresholded clusters)
%      * summary tables with atlas-based annotations
%
% -------------------------------------------------------------------------
% NOTES
% -------------------------------------------------------------------------
% - All statistical inference is non-parametric and permutation-based.
% - TFCE is used as the primary method for spatial inference, with explicit
%   control of family-wise error via the permutation-derived max-statistic.
% - The script is fully non-interactive and suitable for HPC execution.
% - No temporary decoding folders or files are written during permutation
%   decoding.
%
% -------------------------------------------------------------------------
% AUTHOR
% -------------------------------------------------------------------------
% Lukas Van Oudenhove
% Leuven, Belgium
% 2026
%
%
%% ========================================================================
% 0. USER SETTINGS — EDIT THESE
% ========================================================================

% INPUT DIRECTORIES

discoverie_secondlevel_m11_s0_a_set_up_paths_always_run_first
    
% PHENOTYPE FILE WITH GROUP INFO

phenofile = fullfile(BIDSdir, 'phenotype.csv');
phenotype = readtable(phenofile);

% CROSS-VALIDATION SETTINGS

K = 5;

% TDT SETTINGS

con2use = 'con_0003.nii';
conname = DSGN.contrastnames{sscanf(con2use, 'con_%d.nii')};
conname = strrep(conname, ' ', '_');
performance_metric = {'AUC_minus_chance'}; % script only tested with one metric
unbalanced = true;

analysis_mode = 'searchlight';    % 'searchlight' | 'roi' | 'wholebrain'

switch analysis_mode
    case 'searchlight'
        searchlight_radius = 6; % in mm
        mask_file = which('gm_mask_canlab2023_coarse_fmriprep20_0_20.nii');
        mask_obj = fmri_mask_image(mask_file);
    case 'roi'
        roi_list = {
                   which('aINS_L.nii')
                   which('aINS_R.nii')
                   which('Tha_L.nii')
                   which('Tha_R.nii')
                   which('Amyg_L.nii')
                   which('Amyg_R.nii')
                   which('IFG_vlPFC_L.nii')  
                   which('IFG_vlPFC_R.nii')
               }; % Add more ROI masks here - can be created with LaBGAS script
    case 'wholebrain'
        mask_file = which('gm_mask_canlab2023_coarse_fmriprep20_0_20.nii');
    otherwise
        error('analysis mode specified not implemented in TDT')
end

% PERMUTATION SETTINGS

n_perms = 1000;

% OUTPUT DIRECTORIES

tdt_resultsdir = fullfile(resultsdir,'TDT',conname,analysis_mode);  % folder where all TDT results are saved, can work with subfolders if you want to run multiple analyses
    if ~exist(tdt_resultsdir,'dir')
        mkdir(tdt_resultsdir);
    end
    
% ATLAS SETTINGS
atlas = load_atlas('canlab2024');

% THRESHOLDS FOR OUTPUT FMRI_DATA OBJECTS
unc_p = 0.01;
unc_k = 50;
fdr_p = 0.05;
fdr_k = 25;
fwe_p = 0.05;
fwe_k = 25;

    
%% ========================================================================
% 1. BUILD CV DESIGN
% ========================================================================

HC = phenotype(phenotype.patient == 0,:);
HC.id = (1:height(HC))';

IBS = phenotype(phenotype.patient == 1,:);
IBS.id = (1:height(IBS))';

% BUILD CROSS-VALIDATION DESIGN

folds_IBS = cvpartition(height(IBS), 'KFold', K);
folds_HC = cvpartition(height(HC), 'KFold', K);

cv = struct;
cv.HC = folds_HC.TestSize; % 5 fold cross validation, balanced HCs and patients in each fold
cv.IBS = folds_IBS.TestSize;

cv.rand_HC = randperm(height(HC)); % random indices for HCs
cv.rand_IBS = randperm(height(IBS)); % random indices for IBS

% Assign subsets
startIdx = 1;
for i = 1:length(cv.HC) % no return sampling for HCs
    endIdx = startIdx + cv.HC(i) -1;
    cv.subsets_HC{i} = cv.rand_HC(startIdx:endIdx);
    startIdx = endIdx + 1;
end

startIdx = 1;
for i = 1:length(cv.IBS) % no return sampling for IBS
    endIdx = startIdx + cv.IBS(i) -1;
    cv.subsets_IBS{i} = cv.rand_IBS(startIdx:endIdx);
    startIdx = endIdx + 1;
end

% Assign chunk numbers
for k = 1:length(cv.HC) % filling out the chunks for design matrix
    found = ismember(HC.id, cv.subsets_HC{1,k});
    HC.chunks(found) = k;
    
    found = ismember(IBS.id, cv.subsets_IBS{1,k});
    IBS.chunks(found) = k;
end

table_combined = [HC; IBS];
table_combined.id = (1:height(table_combined))';


%% ========================================================================
% 2. BUILD TDT CONFIGURATION (Mode‑Adaptive) - EDIT THESE
% ========================================================================

cfg = decoding_defaults;
cfg.results.dir = tdt_resultsdir;

switch lower(analysis_mode)

   case 'searchlight'
       cfg.analysis = 'searchlight';
       cfg.searchlight.unit = 'mm';
       cfg.searchlight.radius = searchlight_radius;
       cfg.searchlight.spherical = 1;
       cfg.files.mask = mask_file;
       fprintf('>> RUNNING SEARCHLIGHT\n');

   case 'roi'
       cfg.analysis = 'ROI';
       cfg.files.mask = roi_list;
       fprintf('>> RUNNING ROI DECODING\n');

   case 'wholebrain'
       cfg.analysis = 'wholebrain';
       cfg.files.mask = mask_file;
       fprintf('>> RUNNING WHOLE‑BRAIN DECODING\n');

   otherwise
       error('Unknown analysis mode.');
end

% Load subject images
cfg.files.name = cell(height(table_combined),1);
for i = 1:height(table_combined)
   cfg.files.name{i} = sprintf(['%s/%s/' con2use], ...
       datadir, table_combined.PPID{i});
end

cfg.files.chunk = table_combined.chunks;

table_combined.labels(table_combined.patient==1) = 1;
table_combined.labels(table_combined.patient==-1) = -1;
cfg.files.label = table_combined.labels;

cfg.results.output = performance_metric;
cfg.verbose = 0;

% z-scoring training data per fold, apply to test data in same fold
% use this for real decoding in cases with different sites, scanner, etc
cfg.decoding_method = 'classification_kernel';
cfg.scale.method    = 'z';
cfg.scale.estimation = 'separate';

cfg.plot_selected_voxels = 0;

cfg.design = make_design_cv(cfg);
if unbalanced
    cfg.design.unbalanced_data = 'ok';
end


%% ========================================================================
% 3. RUN TRUE DECODING
% ========================================================================

fprintf('\n=== RUNNING TRUE DECODING ===\n');
results = decoding(cfg);

% Extract vector format of real results
outname = cfg.results.output{1};
real_vec = results.(outname).output(:);


%% ========================================================================
% 4. GENERATE PERMUTATION DESIGNS
% ========================================================================

cfgp = rmfield(cfg,'design');
cfgp.design.function.name = 'make_design_cv';
% z-scoring once for all data is acceptable for permutation data as there
% are not scanner/site differences here. 'separate' not viable in terms of
% runtime
cfgp.decoding.method = 'classification';
cfgp.scale.method = 'z';
cfgp.scale.estimation = 'all';

designs = make_design_permutation(cfgp, n_perms, 0);


%% ========================================================================
% 5. FAST PERMUTATION DECODING (NO FOLDERS, IN‑MEMORY)
% ========================================================================

fprintf('\n=== RUNNING %d PERMUTATIONS WITHOUT WRITING FILES ===\n', n_perms);

% Start pool if needed
LaBGAScore_smart_parallel_pool_setup

Cdesigns = parallel.pool.Constant(designs);

% Build template cfg
cfg_template = decoding_defaults;

    % --- Copy ONLY harmless, structural fields ---
    cfg_template.analysis   = cfg.analysis;
    cfg_template.files      = cfg.files;
    cfg_template.searchlight = cfg.searchlight;

    % --- FORCE linear decoding ---
    cfg_template.decoding.method = 'classification';
    cfg_template.decoding.train.classification.model_parameters = '-s 2 -c 1 -q';

    % --- FORCE fixed scaling ---
    cfg_template.scale.method    = 'z';
    cfg_template.scale.estimation = 'all';

    % --- Outputs ---
    cfg_template.results.output = cfg.results.output;
    cfg_template.results.write = 0;
    cfg_template.verbose       = 2;
    cfg_template.plot_selected_voxels = 0;

    cfg_template.design = [];

% Preallocate null matrix
V = numel(real_vec);
all_perm_results = zeros(V, n_perms, 'single');

% Progress
tracker = ProgressTracker(n_perms);
q = parallel.pool.DataQueue;
afterEach(q, @(~) tracker.update());

fprintf('Running %d permutations...\n', n_perms);
disp(cfg_template.searchlight)

% PARFOR: decode permutations in memory only
parfor p = 1:n_perms
   cfgPi = cfg_template;
   cfgPi.design = Cdesigns.Value{p};
   cfgPi.design.unbalanced_data = 'ok';
   
   r = decoding(cfgPi);

   all_perm_results(:,p) = single(r.(outname).output(:));

   send(q,1);
end

save(fullfile(tdt_resultsdir,'combined_permutations.mat'), ...
    'all_perm_results','outname','-v7.3');

fprintf('\nCombined null distribution saved.\n');


%% ========================================================================
% 6. VOXELWISE P-VALUES
% ========================================================================

fprintf('\n=== COMPUTING VOXELWISE P-VALUES ===\n');

[V,P] = size(all_perm_results);

exceed = sum(all_perm_results >= real_vec, 2);
p_unc = (exceed + 1) ./ (P + 1);

switch analysis_mode
    
    case 'wholebrain'
        
        results.(performance_metric{1}).p_perm = p_unc;
        save(fullfile(tdt_resultsdir,['res_' performance_metric{1} '.mat']),'results','-append');
        
    case 'roi'
        
        results.(performance_metric{1}).p_perm_unc = p_unc;
        
        % Apply Storey's FDR
        [~, q_fdr, aprioriprob] = mafdr(p_unc);

        % If aprioriprob > 0.99, fallback to BenjaminiHochberg
        if aprioriprob > 0.99
            p_fdr = mafdr(p_unc, 'BHFDR', true);
            results.(performance_metric{1}).p_perm_fdr = p_fdr;
        else
            % Enforce constraint q >= p (as in SAS proc multtest)
            for j = 1:length(q_fdr)
                if q_fdr(j) < p_unc(j)
                    q_fdr(j) = p_unc(j);
                end
            end
            results.(performance_metric{1}).q_perm_fdr = q_fdr;
        end
        
        save(fullfile(tdt_resultsdir,['res_' performance_metric{1} '.mat']),'results','-append');
        
        
    case 'searchlight'

        results.(performance_metric{1}).p_perm = p_unc;
        
        % Save uncorrected whole-brain p-map
        nii = load_nii(mask_file);
        mask = nii.img > 0;
        mask_idx = find(mask);

        % real_vec must be the real decoding statistics vector
        % used to generate p_unc (same length)
        if length(real_vec) ~= length(p_unc)
            error('real_vec and p_unc length mismatch');
        end

        % Identify valid decoding voxels (NaNs removed by TDT)
        valid_idx = ~isnan(real_vec);

        % Initialize p-map
        pmap = nan(size(nii.img), 'single');

        % Assign only to valid voxels
        pmap(mask_idx(valid_idx)) = p_unc;

        nii.img = pmap;
        save_nii(nii, fullfile(tdt_resultsdir,'p_uncorrected.nii'));
        
        [AUC_stat_obj, AUC_fmri_data, AUC_region_obj, AUC_region_table] = thresholded_fmri_data_from_pval_nii(fullfile(tdt_resultsdir,'p_uncorrected.nii'), results.(performance_metric{1}).output, ...
            mask_obj, atlas, unc_p, ['right-tailed uncorrected p-values based on ' num2str(n_perms) ' permutations'], 'AUC', 'unc', unc_k);


        %% ========================================================================
        % 7. WHOLE-BRAIN FDR CORRECTION
        % ========================================================================

        fprintf('\n=== APPLYING FDR CORRECTION ===\n');
        
        fmap = nan(size(nii.img));

        % Apply Storey's FDR
        [~, q_fdr, aprioriprob] = mafdr(p_unc);

        % If aprioriprob > 0.99, fallback to BenjaminiHochberg
        if aprioriprob > 0.99
            p_fdr = mafdr(p_unc, 'BHFDR', true);
            results.(performance_metric{1}).p_perm_fdr = p_fdr;
            fmap(mask_idx(valid_idx)) = p_fdr;
            nii.img = fmap;
            save_nii(nii, fullfile(tdt_resultsdir,'p_FDR.nii'));
            [AUC_stat_obj_fdr, AUC_fmri_data_fdr, AUC_region_obj_fdr, AUC_region_table_fdr] = thresholded_fmri_data_from_pval_nii(fullfile(tdt_resultsdir,'p_FDR.nii'), results.(performance_metric{1}).output, ...
                mask_obj, atlas, fdr_p, ['right-tailed FDR-corrected p-values based on ' num2str(n_perms) ' permutations'], 'AUC', 'fdr', fdr_k);
        else
            % Enforce constraint q >= p (as in SAS proc multtest)
            for j = 1:length(q_fdr)
                if q_fdr(j) < p_unc(j)
                    q_fdr(j) = p_unc(j);
                end
            end
            results.(performance_metric{1}).q_perm_fdr = q_fdr;
            fmap(mask_idx(valid_idx)) = q_fdr;
            nii.img = fmap;
            save_nii(nii, fullfile(tdt_resultsdir,'q_FDR.nii'));
            [AUC_stat_obj_fdr, AUC_fmri_data_fdr, AUC_region_obj_fdr, AUC_region_table_fdr] = thresholded_fmri_data_from_pval_nii(fullfile(tdt_resultsdir,'q_FDR.nii'), results.(performance_metric{1}).output, ...
                mask_obj, atlas, fdr_p, ['right-tailed FDR-corrected q-values based on ' num2str(n_perms) ' permutations'], 'AUC', 'fdr', fdr_k);
        end

        %% ========================================================================
        % 8. TFCE AND TFCE P-VALUES
        % ========================================================================

        fprintf('\n=== RUNNING pTFCE on t-maps ===\n');

        % -------------------- TFCE parameters --------------------
        voxel_size = double(nii.hdr.dime.pixdim(2:4)); % [vx vy vz]
        H   = 2;         % Height exponent (Smith & Nichols default)
        E   = 0.5;       % Extent exponent
        conn = 26;        % 3D connectivity

        % ========================================================================
        % 0. COMPUTE REAL AND PERMUTATION T-MAPS + EPSILON
        % ========================================================================

        eps_val = 1e-6;     % prevents zero-plateaus that trigger warnings

        % Real t-map
        real_t_vec = (real_vec - mean(all_perm_results,2)) ./ std(all_perm_results,0,2);
        real_t_vec = real_t_vec + eps_val;
        results.(performance_metric{1}).real_t = real_t_vec;

        % Permutation t-maps matrix: [V x P]
        perm_t_vec = (all_perm_results - mean(all_perm_results,2)) ./ std(all_perm_results,0,2);
        perm_t_vec = perm_t_vec + eps_val;
        results.(performance_metric{1}).perm_t = perm_t_vec;

        % ========================================================================
        % 1. REAL TFCE MAP
        % ========================================================================

        tmpVol = zeros(size(nii.img), 'single');
        tmpVol(mask_idx(valid_idx)) = single(real_t_vec);

        try
           real_TFCE_img = pTFCE(tmpVol, voxel_size, H, E, conn);
        catch
           real_TFCE_img = pTFCE(tmpVol, voxel_size, H, E);
        end
        
        real_TFCE_img(isnan(real_TFCE_img)) = 0;
        real_TFCE_img(real_TFCE_img < 0) = 0;

        real_TFCE_vec = real_TFCE_img(mask_idx(valid_idx));
        results.(performance_metric{1}).real_TFCE = real_TFCE_vec;

        % Save TFCE-enhanced real map
        tfmap = nan(size(nii.img), 'single');
        tfmap(mask_idx(valid_idx)) = real_TFCE_vec;
        nii.img = tfmap;
        save_nii(nii, fullfile(tdt_resultsdir,'TFCE_real.nii'));

        % ========================================================================
        % 2. PARALLEL TFCE NULL DISTRIBUTION (t-maps → TFCE)
        % ========================================================================

        fprintf('Computing TFCE null distribution (parallel)...\n');

        TFCE_null     = zeros(P,1,'single');     % global max TFCE
        TFCE_perm_all = zeros(V,P,'single');     % voxel-wise TFCE
        
        % Progress
        tracker = ProgressTracker(P);
        q_tfce = parallel.pool.DataQueue;
        afterEach(q_tfce, @(~) tracker.update());

        parfor p = 1:P

           tmpVol = zeros(size(nii.img), 'single');
           tmpVol(mask_idx(valid_idx)) = single(perm_t_vec(:,p));

           % Compute TFCE on permutation t-map
           try
               TFCE_img = pTFCE(tmpVol, voxel_size, H, E, conn);
           catch
               TFCE_img = pTFCE(tmpVol, voxel_size, H, E);
           end
           
           TFCE_img(isnan(TFCE_img)) = 0;

           TFCE_null(p) = max(TFCE_img(:));         % global max-stat
           TFCE_perm_all(:,p) = TFCE_img(mask_idx(valid_idx)); % voxel-wise TFCE
           
           send(q_tfce,1);
           
        end
        
        TFCE_perm_all(TFCE_perm_all < 0) = 0;

        save(fullfile(tdt_resultsdir,'TFCE_null.mat'),'TFCE_null','TFCE_perm_all');

        % ========================================================================
        % 3. GLOBAL FWE-CORRECTED TFCE p-VALUE (max-statistic)
        % ========================================================================

        TFCE_real_max = max(real_TFCE_vec);

        p_TFCE_global = (sum(TFCE_null >= TFCE_real_max) + 1) / (P + 1);
        
        results.(performance_metric{1}).p_TFCE_global = p_TFCE_global;

        % ========================================================================
        % 4. VOXEL-WISE TFCE p-MAP
        % ========================================================================

        fprintf('Computing voxel-wise TFCE p-values...\n');

        p_TFCE_voxelwise = (sum(TFCE_perm_all >= real_TFCE_vec, 2) + 1) / (P + 1);
        
        results.(performance_metric{1}).p_TFCE_unc = p_TFCE_voxelwise;

        pmap_vox = nan(size(nii.img), 'single');
        pmap_vox(mask_idx(valid_idx)) = p_TFCE_voxelwise;

        nii.img = pmap_vox;
        save_nii(nii, fullfile(tdt_resultsdir,'p_TFCE_voxelwise.nii'));
        
        [TFCE_stat_obj, TFCE_fmri_data, TFCE_region_obj, TFCE_region_table] = thresholded_fmri_data_from_pval_nii(fullfile(tdt_resultsdir,'p_TFCE_voxelwise.nii'), results.(performance_metric{1}).real_TFCE, ...
            mask_obj, atlas, unc_p, ['right-tailed uncorrected TFCE p-values based on ' num2str(n_perms) ' permutations'], 'TFCE', 'unc', unc_k);

        % ========================================================================
        % 5. FDR-CORRECTED TFCE p-MAP
        % ========================================================================

        fprintf('Applying FDR correction to TFCE p-values...\n');
        
        fmap = nan(size(nii.img));

        % Apply Storey's FDR
        [~, q_TFCE_FDR, aprioriprob] = mafdr(p_TFCE_voxelwise);

        % If aprioriprob > 0.99, fallback to Benjamini-Hochberg
        if aprioriprob > 0.99
            p_TFCE_FDR = mafdr(p_TFCE_voxelwise, 'BHFDR', true);
            results.(performance_metric{1}).p_TFCE_FDR = p_TFCE_FDR;
            pmap_FDR = nan(size(nii.img), 'single');
            pmap_FDR(mask_idx(valid_idx)) = p_TFCE_FDR;
            nii.img = pmap_FDR;
            save_nii(nii, fullfile(tdt_resultsdir,'p_TFCE_FDR_voxelwise.nii'));
            
            [TFCE_stat_obj_fdr, TFCE_fmri_data_fdr, TFCE_region_obj_fdr, TFCE_region_table_fdr] = thresholded_fmri_data_from_pval_nii(fullfile(tdt_resultsdir,'p_TFCE_FDR_voxelwise.nii'), results.(performance_metric{1}).real_TFCE, ...
                mask_obj, atlas, fdr_p, ['right-tailed fdr-corrected TFCE p-values based on ' num2str(n_perms) ' permutations'], 'TFCE', 'fdr', fdr_k);

        else
            % Enforce constraint q >= p (as in SAS proc multtest)
            for j = 1:length(q_TFCE_FDR)
                if q_TFCE_FDR(j) < p_TFCE_voxelwise(j)
                    q_TFCE_FDR(j) = p_TFCE_voxelwise(j);
                end
            end
            results.(performance_metric{1}).q_TFCE_FDR = q_TFCE_FDR;
            qmap_FDR = nan(size(nii.img), 'single');
            qmap_FDR(mask_idx(valid_idx)) = q_TFCE_FDR;
            nii.img = qmap_FDR;
            save_nii(nii, fullfile(tdt_resultsdir,'q_TFCE_FDR_voxelwise.nii'));
            
            [TFCE_stat_obj_fdr, TFCE_fmri_data_fdr, TFCE_region_obj_fdr, TFCE_region_table_fdr] = thresholded_fmri_data_from_pval_nii(fullfile(tdt_resultsdir,'q_TFCE_FDR_voxelwise.nii'), results.(performance_metric{1}).real_TFCE, ...
                mask_obj, atlas, fdr_p, ['right-tailed fdr-corrected TFCE q-values based on ' num2str(n_perms) ' permutations'], 'TFCE', 'fdr', fdr_k);

        end
        
        %% ========================================================================
        % 6. VOXEL-WISE TFCE FWE-CORRECTED P-MAP (MAX-STATISTIC)
        % ========================================================================

        fprintf('Computing voxel-wise TFCE FWE-corrected p-values (max-statistic)...\n');

        P = numel(TFCE_null); % number of permutations

        % -------------------------------------------------------------------------
        % Compute voxel-wise TFCE FWE p-values using max TFCE null
        % -------------------------------------------------------------------------
        p_TFCE_FWE_voxelwise = ( ...
           sum(TFCE_null(:)' >= real_TFCE_vec, 2) + 1 ...
           ) ./ (P + 1);
        results.(performance_metric{1}).p_TFCE_FWE = p_TFCE_FWE_voxelwise;
        
        save(fullfile(tdt_resultsdir,['res_' performance_metric{1} '.mat']),'results','-append');

        % -------------------------------------------------------------------------
        % Map vector -> volume
        % -------------------------------------------------------------------------
        pmap_FWE = nan(size(nii.img), 'single');
        pmap_FWE(mask_idx(valid_idx)) = p_TFCE_FWE_voxelwise;

        nii.img = pmap_FWE;
        save_nii(nii, fullfile(tdt_resultsdir,'p_TFCE_FWE_voxelwise.nii'));
        
        [TFCE_stat_obj_fwe, TFCE_fmri_data_fwe, TFCE_region_obj_fwe, TFCE_region_table_fwe] = thresholded_fmri_data_from_pval_nii(fullfile(tdt_resultsdir,'p_TFCE_FWE_voxelwise.nii'), results.(performance_metric{1}).real_TFCE, ...
            mask_obj, atlas, fwe_p, ['right-tailed fwe-corrected TFCE p-values (max statistic) based on ' num2str(n_perms) ' permutations'], 'TFCE', 'fwe', fwe_k);
        
        fprintf('Voxel-wise TFCE FWE-corrected map saved.\n');

        fprintf('\nTFCE completed: FWE, voxel-wise, and FDR TFCE maps saved.\n');

        
        save(fullfile(tdt_resultsdir,'stat_objs.mat'),'AUC_*','TFCE_*');

        
end

%% HELPER FUNCTION

function [stat_img, fmri_dat, region_obj, region_table] = thresholded_fmri_data_from_pval_nii (fullpath_to_pval_nii, stat_vec, mask_obj, atlas_obj, p, p_type, stat_type, thr_type, k)

        stat_img = statistic_image(fullpath_to_pval_nii,'type','p');
        stat_img = apply_mask(stat_img,mask_obj);
        stat_img.sig = logical(stat_img.p < p); % manual thresholding at 0.05, statistic_image.threshold() breaks for some reason
        stat_img.p_type = p_type;
        
        if sum(stat_img.sig) == 0
            fprintf('\n No suprathreshold voxels after p threshold\n');
            fmri_dat = [];
            region_obj = [];
            region_table = [];
            
        else
        
            fmri_dat = fmri_data_st(stat_img);
            fmri_dat.dat = stat_vec;
            fmri_dat.dat_descrip = [stat_type ' thresholded at p_' thr_type ' < ' num2str(p)];
            path = fileparts(fullpath_to_pval_nii); 
            fmri_dat = fmri_dat.threshold([min(fmri_dat.dat(stat_img.sig))-0.001 max(fmri_dat.dat(stat_img.sig))+0.001],'raw-between','k',k);
            fmri_dat.image_names = [stat_type '_' strrep(num2str(p),'.','_') '_' thr_type '_k_' num2str(k) '.nii'];
            fmri_dat.fullpath = fullfile(path,[stat_type '_' strrep(num2str(p),'.','_') '_' thr_type '_k_' num2str(k) '.nii']);
            
            if isempty(fmri_dat)
                fprintf('\n No suprathreshold voxels after k threshold\n');
                region_obj = [];
                region_table = [];
            
            else

                region_obj = region(fmri_dat);
                region_obj = region_obj.autolabel_regions_using_atlas(atlas_obj);

                [~,~,region_table] = table(region_obj);
                average_region = zeros(height(region_table),1);
                k_region = zeros(height(region_table),1);
                for r = 1:size(region_obj,2)
                    average_region(r,1) = mean(region_obj(1,r).val);
                    k_region(r,1) = region_obj(1,r).numVox;
                end
                region_table.minP = [];
                region_table.(['mean_' stat_type]) = average_region;
                region_table.k = k_region;
                
            end
            
        end
        
end