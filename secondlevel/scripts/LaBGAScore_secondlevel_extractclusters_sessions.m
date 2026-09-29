%% LaBGAScore_secondlevel_extractclusters_sessions.m
%
%
% *USAGE*
%
% The voxel-wise counterpart of LaBGAScore_secondlevel_extractparcels_sessions.
% Takes the CLUSTERS that survived a BETWEEN-session contrast in the voxel-wise
% GLM, extracts each subject's mean value within each cluster for the
% WITHIN-session contrasts, and builds a long-format table (one row per subject
% per session) ready for mixed models, plotting or export.
%
% Run it after BOTH prep_3a_run_second_level_regression_and_save.m AND
% c2a_second_level_regression.m: the clusters are written by c2a, not prep_3a.
% c2a APPENDS its region objects into the same
% regression_stats_and_maps_<...>.mat that prep_3a wrote, so both the fit and
% the clusters come from one file.
%
%
% *WHY THIS IS NOT THE SAME OPERATION AS extractparcels*
%
% A PARCEL is a column of parcelwise_stats_results.datmatrix, so the parcel
% script extracts by index - the per-subject values already exist.
%
% A CLUSTER is a set of voxels, and the region objects c2a saves carry geometry
% and statistic ONLY:
%
%   XYZ       3 x numVox     voxel coordinates
%   Z         1 x numVox     the statistic per voxel
%   val       numVox x 1
%   dat       0 x 0          EMPTY
%   all_data  0 x 0          EMPTY
%
% c2a builds them from thresholded statistic images, which have no subject
% dimension, so there is nothing per-subject to look up. The values must be
% extracted SPATIALLY: each subject's within-session contrast image averaged
% over the cluster's voxels. This script therefore needs a second input,
% contrast_data_objects.mat, and is slower than the parcel version.
%
%
% *OPTIONS*
%
% * cluster_type                       'fdr' | 'tfce_fdr' | 'tfce_fwe' | 'unc'
% * nr_sess                            number of within-subject sessions
% * names_sess                         session labels, one per session
% * within_session_contrast_idx        indices from DAT.contrastnames for the within-session contrasts
% * between_session_contrast_idx       index from DAT.contrastnames for the between-session contrast
% * mygroupnamefield                   'conditions' or 'contrasts', must match the prep_3a script
% * results_suffix                     results_suffix of the VOXEL-WISE prep_3a run
% * dosave                             write the table to resultsdir as .csv, default true
%
%
% *OUTPUT*
%
% * results_table_cluster   one row per subject per session, one column per
%                           cluster, plus PPID, group, session and session_name.
%                           Also written to resultsdir as .csv.
%
%
% *THREE TRAPS, ALL OF WHICH FAIL SILENTLY*
%
% Recorded because each returns wrong output rather than an error:
%
% 1. @fmri_data/extract_roi_averages REJECTS A REGION ARRAY ("unknown mask
%    input type") - it takes a mask image. @region/extract_data is the method
%    for this. It works from mm coordinates, so the data object need not share
%    the statistic image's space, and it fills r(i).dat with the region average
%    per subject and r(i).all_data with the voxelwise data.
%
% 2. region(img,'unique_mask_values') RETURNS ONE EMPTY REGION on a continuous
%    statistic map. That flag groups voxels by identical intensity and assumes
%    an integer-coded atlas. Use the DEFAULT, 'contiguous_regions', which is
%    what c2a itself uses (a bare region(tfce_dat_thr_corr)). Measured side by
%    side on one FDR-thresholded TFCE image: unique_mask_values gave 1 empty
%    region with numVox = 0; the default gave 193.
%
% 3. region_objs_tfce_corr HOLDS ONLY ONE CORRECTION - whichever
%    tfce_correction names. Asking for TFCE-FDR while it holds FWE and trusting
%    it returns FWE clusters under an FDR label. Measured on one model: 116 FWE
%    clusters where the FDR answer was 193. This script reads tfce_correction
%    and, on a mismatch, rebuilds from the corresponding thresholded statistic
%    image (tfce_stat_imgs_thr_fdr / _fwe) and relabels via the atlas so the
%    names stay anatomical rather than Region001.
%
% A related trap when writing any voxel extraction by hand: volInfo.wh_inmask
% is NOT a 1:1 index into obj.dat. CANlab drops voxels from .dat without
% shrinking the index (235807 rows against 235820 entries in one measured case),
% so indexing .dat by position in wh_inmask returns wrong values with no error.
% @region/extract_data does its own coordinate matching and avoids this.
%
%
% *DEPENDENCIES*
%
% CanlabCore (@region/extract_data, region, autolabel_regions_using_atlas,
% printhdr), and the .mat result files produced by prep_3a and c2a. Requires
% myscaling_glm already in the workspace (normally via a2_set_default_options.m /
% a_set_up_paths_always_run_first).
%
% -------------------------------------------------------------------------
%
% author: Lukas Van Oudenhove
%
% date:   KU Leuven, September 2026
%
% -------------------------------------------------------------------------
%
% LaBGAScore_secondlevel_extractclusters_sessions.m          v1.0
%
% last modified: 2026/09/29


%% ========================================================================
% 0. USER SETTINGS — EDIT THESE
% =========================================================================

clear cluster_names col_idx

cluster_type = 'tfce_fwe';             % 'fdr' | 'tfce_fdr' | 'tfce_fwe' | 'unc'
dosave       = true;

% INPUT DIRECTORIES

LaBGAScore_prep_s0_define_directories;
a_set_up_paths_always_run_first;
load(fullfile(resultsdir,'image_names_and_setup.mat'));

% SESSION INFO

nr_sess                      = 2;
names_sess                   = {'ses-01','ses-02'};
within_session_contrast_idx  = [1,2];  % indices from DAT.contrastnames
between_session_contrast_idx = 5;

% MUST MATCH THE VOXEL-WISE prep_3a RUN THAT c2a REPORTED ON

mygroupnamefield = 'contrasts';
results_suffix   = '';                 % e.g. 'vox_gm'

switch myscaling_glm

    case 'raw'
        fprintf('\nContrast calculated on raw (unscaled) condition images used in second-level GLM\n\n');
        scaling_string = 'no_scaling';          dataobj_name = 'DATA_OBJ_CON';

    case 'scaled'
        fprintf('\nContrast calculated on z-scored condition images used in second-level GLM\n\n');
        scaling_string = 'scaling_z_score_conditions';  dataobj_name = 'DATA_OBJ_CONsc';

    case 'scaled_contrasts'
        fprintf('\nl2norm scaled contrast images used in second-level GLM\n\n');
        scaling_string = 'scaling_l2norm_contrasts';    dataobj_name = 'DATA_OBJ_CONscc';

    otherwise
        error(['Invalid option "%s" in myscaling_glm (a2_set_default_options): ' ...
               'choose "raw", "scaled" or "scaled_contrasts".'], myscaling_glm);

end

fprintf('scaling %s -> extracting from %s\n\n', myscaling_glm, dataobj_name);

if numel(within_session_contrast_idx) ~= nr_sess
    error('within_session_contrast_idx has %d entries but nr_sess is %d.', ...
        numel(within_session_contrast_idx), nr_sess);
end
if numel(names_sess) ~= nr_sess
    error('names_sess has %d entries but nr_sess is %d.', numel(names_sess), nr_sess);
end


%% ========================================================================
% 1. LOAD THE CLUSTERS
% =========================================================================

fprintf('\n\n'); printhdr('LOADING CLUSTERS'); fprintf('\n\n');

f_reg = fullfile(resultsdir, ['regression_stats_and_maps_', mygroupnamefield, '_', ...
                              scaling_string, '_', results_suffix, '.mat']);
if ~exist(f_reg,'file')
    error(['No voxel-wise results at %s.\nRun prep_3a and then c2a - the clusters ' ...
           'are written by c2a, not prep_3a.'], f_reg);
end
fprintf('reading clusters from %s\n', f_reg);

V    = whos('-file', f_reg);
have = @(n) any(strcmp({V.name}, n));

switch lower(cluster_type)

    case 'unc'
        if ~have('region_objs_unc'), error('region_objs_unc not in %s; has c2a been run?', f_reg); end
        Q = load(f_reg,'region_objs_unc');
        regs = Q.region_objs_unc{between_session_contrast_idx};
        src  = 'region_objs_unc';

    case 'fdr'
        if ~have('region_objs_fdr'), error('region_objs_fdr not in %s; has c2a been run?', f_reg); end
        Q = load(f_reg,'region_objs_fdr');
        regs = Q.region_objs_fdr{between_session_contrast_idx};
        src  = 'region_objs_fdr';

    case {'tfce_fdr','tfce_fwe'}
        want = extractAfter(lower(cluster_type), 'tfce_');
        if ~have('tfce_correction')
            error(['No TFCE results in %s. Either c2a ran without TFCE, or dotfce ' ...
                   'was false in prep_3a.'], f_reg);
        end
        Q      = load(f_reg,'tfce_correction');
        stored = lower(strtrim(Q.tfce_correction));
        fprintf('tfce_correction stored in the results = "%s"; requested "%s"\n', stored, want);

        if strcmp(stored, want)
            Q2   = load(f_reg,'region_objs_tfce_corr');
            regs = Q2.region_objs_tfce_corr{between_session_contrast_idx};
            src  = sprintf('region_objs_tfce_corr (tfce_correction = %s)', stored);
        else
            % Trap 3: the stored regions are the OTHER correction. Rebuild
            % rather than silently returning the wrong set.
            imgs_var = sprintf('tfce_stat_imgs_thr_%s', want);
            if ~have(imgs_var)
                error(['Requested TFCE %s, but the stored region objects are %s and %s ' ...
                       'is not in the results file. Re-run c2a with tfce_correction = ' ...
                       '''%s'' to get them.'], upper(want), upper(stored), imgs_var, want);
            end
            fprintf('stored regions are %s, so rebuilding %s clusters from %s\n', ...
                    upper(stored), upper(want), imgs_var);
            Q2  = load(f_reg, imgs_var);
            img = Q2.(imgs_var){between_session_contrast_idx};
            if isempty(img)
                regs = [];
            else
                % Trap 2: the DEFAULT grouping is 'contiguous_regions', which is
                % what c2a uses. Do not pass 'unique_mask_values'.
                regs = region(img);
                if ismember('autolabel_regions_using_atlas', methods('region'))
                    try
                        regs = autolabel_regions_using_atlas(regs);
                    catch ME
                        warning('%s', ['Could not autolabel the rebuilt regions (' ...
                                ME.message '); names stay generic.']);
                    end
                end
            end
            src = sprintf('%s (rebuilt, contiguous regions)', imgs_var);
        end

    otherwise
        error('Invalid cluster_type "%s": choose fdr, tfce_fdr, tfce_fwe or unc.', cluster_type);

end

if iscell(regs), regs = regs{1}; end
n_clu = numel(regs);

fprintf('\nsource: %s\n', src);
fprintf('contrast %d ("%s"): %d cluster(s)\n\n', between_session_contrast_idx, ...
        DAT.contrastnames{between_session_contrast_idx}, n_clu);

results_table_cluster = table();


%% ========================================================================
% 2. EXTRACT PER-SUBJECT CLUSTER MEANS
% =========================================================================

if n_clu == 0

    fprintf(['Nothing survives %s for this contrast, so there is no table to build.\n' ...
             'That is a result, not a failure.\n'], cluster_type);

else

    fprintf('\n\n'); printhdr('EXTRACTING PER-SUBJECT CLUSTER MEANS'); fprintf('\n\n');

    f_dat = fullfile(resultsdir,'contrast_data_objects.mat');
    if ~exist(f_dat,'file')
        error('No %s. Run prep_3_calc_univariate_contrast_maps_and_save.m first.', f_dat);
    end
    D    = load(f_dat, dataobj_name);
    DOBJ = D.(dataobj_name);

    % Cluster names, made unique and valid. Several clusters routinely share a
    % shorttitle - one anatomical label split into disjoint blobs - and unlike
    % parcels those are DISTINCT clusters that must all be kept, so they are
    % suffixed rather than deduplicated.
    raw = {regs.shorttitle};
    raw(cellfun(@isempty, raw)) = {'unlabeled'};
    cluster_names = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(raw));
    n_dup = numel(raw) - numel(unique(raw));
    if n_dup > 0
        fprintf('%d cluster(s) share a label with another; names suffixed to keep them distinct\n', n_dup);
    end

    tabs = cell(1,nr_sess);
    for s = 1:nr_sess

        obj = DOBJ{within_session_contrast_idx(1,s)};
        fprintf('session %d (contrast %d, "%s"): %d images\n', s, ...
                within_session_contrast_idx(1,s), ...
                DAT.contrastnames{within_session_contrast_idx(1,s)}, size(obj.dat,2));

        % Trap 1: @region/extract_data, NOT @fmri_data/extract_roi_averages.
        ex = extract_data(regs, obj);

        vals = nan(size(obj.dat,2), n_clu);
        for k = 1:n_clu
            v = ex(k).dat;
            if isempty(v)
                error('Cluster %d ("%s") extracted no data.', k, cluster_names{k});
            end
            vals(:,k) = v(:);
        end

        T = array2table(vals, 'VariableNames', cluster_names);
        T.PPID         = DAT.BEHAVIOR.behavioral_data_table.participant_id;
        T.group        = DAT.BEHAVIOR.behavioral_data_table.group;
        T.session      = s .* ones(height(T),1);
        T.session_name = repmat(string(names_sess{s}), height(T), 1);
        tabs{s} = T;

    end

    results_table_cluster = vertcat(tabs{:});
    fprintf('\ncluster table: %d rows x %d columns\n', ...
            height(results_table_cluster), width(results_table_cluster));

end


%% ========================================================================
% 3. SAVE
% =========================================================================

if dosave && ~isempty(results_table_cluster)

    fprintf('\n\n'); printhdr('WRITING TABLE'); fprintf('\n\n');

    fout = fullfile(resultsdir, ['extracted_clusters_', lower(cluster_type), '_', ...
                                 scaling_string, '_', results_suffix, '.csv']);
    writetable(results_table_cluster, fout);
    fprintf('wrote %s\n', fout);

end

fprintf('\n\n'); printhdr('DONE'); fprintf('\n\n');
