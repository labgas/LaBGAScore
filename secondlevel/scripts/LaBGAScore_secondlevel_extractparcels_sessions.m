%% LaBGAScore_secondlevel_extractparcels_sessions.m
%
%
% *USAGE*
%
% This script extracts parcel- and roi-wise WITHIN-session contrast values for
% the parcels that were significant in a BETWEEN-session contrast, and builds
% long-format tables (one row per subject per session) that can be used
% directly for mixed models, plotting, or export.
%
% Run it after the prep_3a script(s) whose results it reads. Note that the
% parcel-wise and roi-wise results are usually written by DIFFERENT prep_3a
% variants and therefore carry DIFFERENT results_suffix values - see
% results_suffix_parcel and results_suffix_roi below.
%
%
% *OPTIONS*
%
% * do_parcel                          run parcel-wise extraction, default true
% * parcel_type                        threshold used to select parcels, 'unc' | 'fdr' | 'Bayes'
% * do_roi                             run roi-wise extraction, default true
% * dosave                             write the tables to resultsdir as .csv, default true
% * nr_sess                            number of within-subject sessions
% * names_sess                         session labels, one per session
% * within_session_contrast_idx        indices from DAT.contrastnames for the within-session contrasts
% * between_session_contrast_idx       index from DAT.contrastnames for the between-session contrast
% * mygroupnamefield                   'conditions' or 'contrasts', must match the prep_3a script
% * results_suffix_parcel              results_suffix of the PARCELWISE prep_3a run
% * results_suffix_roi                 results_suffix of the ROI prep_3a run
%
%
% *OUTPUT*
%
% * results_table_parcel   one row per subject per session, one column per
%                          surviving parcel, plus PPID, group, session and
%                          session_name
% * results_table_roi      the same for the a priori ROIs
%
% Both are also written to resultsdir as .csv.
%
%
% *NOTES*
%
% v2.0 fixes four defects found by running v1.1 against real data. They are
% recorded here because none of them announces itself - three fail silently:
%
% 1. ONE results_suffix CANNOT SERVE BOTH BRANCHES. The parcelwise and roi
%    results are normally written by different prep_3a variants with different
%    suffixes (e.g. 'parc_gm' and 'vox_gm'), so a single suffix finds at most
%    one of the two files. The other branch then printed "no saved results" and
%    returned.
%
% 2. return INSIDE A BRANCH EXITS THE WHOLE SCRIPT. This is a script, not a
%    function, so a missing parcelwise file aborted the roi branch too, even
%    though it is independent and its file was present. The branches now warn
%    and skip.
%
% 3. COLUMN NAMES WERE INDEXED BY REGION COUNT, NOT TABLE WIDTH. v1.1 named
%    columns with a loop over 1:numel(regions) while the table has one column
%    per DISTINCT parcel. Several regions routinely share one parcel label -
%    thresholding splits a parcel into disjoint clusters - so the loop ran past
%    the end of the parcel columns and renamed PPID, group and session to
%    region names before erroring. Silent label/data corruption. Measured on
%    proj_moodbugs_wp2 model_3 with parcel_type = 'Bayes': 115 regions, 21 of
%    them unlabelled and 34 sharing a parcel, giving 60 columns against a loop
%    to 115.
%
% 4. UNMATCHED REGIONS WERE SWEPT AWAY. strcmp against the atlas labels returns
%    all-false for a region whose shorttitle is not a label, and summing that
%    away left fewer parcels than the user believed. The two causes are now
%    separated: a shorttitle of 'No label' is the atlas's marker for a cluster
%    outside every parcel and is DROPPED with a count, which is the correct
%    answer; any other unmatched name means the region set and the atlas
%    disagree (usually atlasname_glm or atlas_granularity not matching the
%    prep_3a run) and is an ERROR.
%
% Also in v2.0: keyword/idx are cleared and pre-allocated, because in a script
% they persist in the base workspace between runs - running 'Bayes' (115
% regions) then 'unc' (1) left 114 stale rows silently selecting the previous
% run's parcels; the tables are written to disk rather than left in the
% workspace only; names_sess is used rather than defined and ignored; and the
% roi branch's printhdr no longer says "PARCELWISE".
%
% The roi branch does NOT add a group column: roi_means_table already carries
% one. The parcel branch does, since its datmatrix does not.
%
%
% *DEPENDENCIES*
%
% Requires myscaling_glm, atlasname_glm and atlas_granularity already in the
% workspace (normally via a2_set_default_options.m /
% a_set_up_paths_always_run_first), and the .mat result files produced by
% prep_3a_run_second_level_regression_and_save.m
%
% -------------------------------------------------------------------------
%
% modified by: Lukas Van Oudenhove
%
% date:   KU Leuven, April 2026
%
% -------------------------------------------------------------------------
%
% LaBGAScore_secondlevel_extractparcels_sessions.m          v2.0
%
% last modified: 2026/09/29


%% ========================================================================
% 0. USER SETTINGS — EDIT THESE
% =========================================================================

clear keyword idx idx_sum col_idx      % see note 5: this is a script, they persist

% CHOOSE ROI AND/OR PARCELS, AND THRESHOLD FOR PARCELS

do_parcel = true;
    parcel_type = 'unc';               % 'unc' | 'fdr' | 'Bayes'
do_roi    = true;
dosave    = true;

% INPUT DIRECTORIES

LaBGAScore_prep_s0_define_directories;
a_set_up_paths_always_run_first;
load(fullfile(resultsdir,'image_names_and_setup.mat'));

% SESSION INFO

nr_sess                      = 2;
names_sess                   = {'ses-01','ses-02'};
within_session_contrast_idx  = [1,2];  % indices from DAT.contrastnames
between_session_contrast_idx = 5;

% SET MANDATORY OPTIONS FROM THE CORRESPONDING PREP_3a SCRIPTS
%
% The parcelwise and roi results normally come from DIFFERENT prep_3a variants
% and therefore carry different suffixes. Set each to the results_suffix of the
% run that produced it.

mygroupnamefield      = 'contrasts';
results_suffix_parcel = '';            % e.g. 'parc_gm'
results_suffix_roi    = '';            % e.g. 'vox_gm'

switch myscaling_glm

    case 'raw'
        fprintf('\nContrast calculated on raw (unscaled) condition images used in second-level GLM\n\n');
        scaling_string = 'no_scaling';

    case 'scaled'
        fprintf('\nContrast calculated on z-scored condition images used in second-level GLM\n\n');
        scaling_string = 'scaling_z_score_conditions';

    case 'scaled_contrasts'
        fprintf('\nl2norm scaled contrast images used in second-level GLM\n\n');
        scaling_string = 'scaling_l2norm_contrasts';

    otherwise
        error(['Invalid option "%s" in myscaling_glm (a2_set_default_options): ' ...
               'choose "raw", "scaled" or "scaled_contrasts".'], myscaling_glm);

end

if numel(within_session_contrast_idx) ~= nr_sess
    error('within_session_contrast_idx has %d entries but nr_sess is %d.', ...
        numel(within_session_contrast_idx), nr_sess);
end
if numel(names_sess) ~= nr_sess
    error('names_sess has %d entries but nr_sess is %d.', numel(names_sess), nr_sess);
end


%% ========================================================================
% 1. PARCEL-WISE EXTRACTION
% =========================================================================

results_table_parcel = table();

if do_parcel

    fprintf('\n\n'); printhdr('LOADING PARCELWISE RESULTS'); fprintf('\n\n');

    savefilenamedata = fullfile(resultsdir, ['parcelwise_stats_and_maps_', ...
        mygroupnamefield, '_', scaling_string, '_', results_suffix_parcel, '.mat']);

    if ~exist(savefilenamedata,'file')
        warning(['No parcelwise results at %s. Skipping the parcel branch.\n' ...
                 'Run the parcelwise prep_3a first, and check results_suffix_parcel.'], savefilenamedata);
        do_parcel = false;
    else
        fprintf('\nLoading parcel-wise results from %s\n\n', savefilenamedata);
        P = load(savefilenamedata);

        atlas            = load_atlas(atlasname_glm);
        atlas_downsample = downsample_parcellation(atlas, ['labels_' num2str(atlas_granularity)]);

        switch parcel_type
            case 'unc',   region_objs = P.region_objs_unc;
            case 'fdr',   region_objs = P.region_objs_fdr;
            case 'Bayes', region_objs = P.region_objs_Bayes;
            otherwise, error('Invalid parcel_type "%s": choose "unc", "fdr" or "Bayes".', parcel_type);
        end

        regions_between_sess = region_objs{1,between_session_contrast_idx}{1,1};
        n_reg = numel(regions_between_sess);

        if n_reg == 0
            warning(['No parcels survive %s thresholding for contrast %d (%s). ' ...
                     'Nothing to extract; skipping the parcel branch.'], ...
                     parcel_type, between_session_contrast_idx, ...
                     DAT.contrastnames{between_session_contrast_idx});
            do_parcel = false;
        else
            fprintf('%d parcel(s) survive %s thresholding for "%s"\n\n', ...
                n_reg, parcel_type, DAT.contrastnames{between_session_contrast_idx});

            % Resolve each region to its atlas column. See notes 3 and 4.
            raw_names  = {regions_between_sess.shorttitle};
            is_nolabel = cellfun(@(x) isempty(x) || strcmpi(strtrim(x),'No label'), raw_names);

            if any(is_nolabel)
                fprintf('dropping %d region(s) with no atlas label (clusters outside every parcel)\n', sum(is_nolabel));
            end

            cand    = raw_names(~is_nolabel);
            keyword = cell(numel(cand),1);
            col_idx = nan(numel(cand),1);
            bad     = {};
            for r = 1:numel(cand)
                keyword{r} = cand{r};
                hits = find(strcmp(atlas_downsample.labels, keyword{r}));
                if isempty(hits)
                    bad{end+1} = keyword{r}; %#ok<SAGROW>
                elseif numel(hits) > 1
                    error('Region "%s" matches %d atlas labels; cannot resolve to one parcel.', ...
                          keyword{r}, numel(hits));
                else
                    col_idx(r) = hits;
                end
            end

            if ~isempty(bad)
                error(['%d region(s) do not match any label in atlas "%s" at granularity %d: %s. ' ...
                       'The region set and the atlas disagree - check that atlasname_glm and ' ...
                       'atlas_granularity match the prep_3a run.'], ...
                       numel(bad), atlasname_glm, atlas_granularity, strjoin(bad, ', '));
            end

            % Several regions can share a parcel: thresholding splits one parcel
            % into disjoint clusters. datmatrix has one column per parcel, so the
            % repeats are redundant - extract each parcel once.
            n_before = numel(col_idx);
            [col_idx, keep] = unique(col_idx, 'stable');
            keyword = keyword(keep);
            if numel(col_idx) < n_before
                fprintf('%d region(s) shared a parcel with another; extracting %d distinct parcel(s)\n', ...
                        n_before - numel(col_idx), numel(col_idx));
            end
            fprintf('extracting %d parcel(s)\n', numel(col_idx));

            results_tables_parcel = cell(1,nr_sess);
            for s = 1:nr_sess
                betas = P.parcelwise_stats_results{1,within_session_contrast_idx(1,s)}.datmatrix;
                if size(betas,2) < max(col_idx)
                    error('datmatrix has %d parcels but a region maps to column %d.', ...
                          size(betas,2), max(col_idx));
                end
                T = array2table(betas(:,col_idx), 'VariableNames', matlab.lang.makeValidName(keyword));
                T.PPID         = DAT.BEHAVIOR.behavioral_data_table.participant_id;
                T.group        = DAT.BEHAVIOR.behavioral_data_table.group;
                T.session      = s .* ones(height(T),1);
                T.session_name = repmat(string(names_sess{s}), height(T), 1);
                results_tables_parcel{s} = T;
            end

            results_table_parcel = vertcat(results_tables_parcel{:});
            fprintf('parcel table: %d rows x %d columns\n', ...
                    height(results_table_parcel), width(results_table_parcel));
        end
    end

end


%% ========================================================================
% 2. ROI-WISE EXTRACTION
% =========================================================================

results_table_roi = table();

if do_roi

    fprintf('\n\n'); printhdr('LOADING ROI RESULTS'); fprintf('\n\n');

    savefilenamedata_roi = fullfile(resultsdir, ['roi_stats_', ...
        mygroupnamefield, '_', scaling_string, '_', results_suffix_roi, '.mat']);

    if ~exist(savefilenamedata_roi,'file')
        warning(['No roi results at %s. Skipping the roi branch.\n' ...
                 'Run the roi prep_3a first, and check results_suffix_roi.'], savefilenamedata_roi);
        do_roi = false;
    else
        fprintf('\nLoading roi-wise results from %s\n\n', savefilenamedata_roi);
        R = load(savefilenamedata_roi);

        results_tables_roi = cell(1,nr_sess);
        for s = 1:nr_sess
            T = R.roi_means_table{1,within_session_contrast_idx(1,s)};
            % roi_means_table already carries group, so it is not added here.
            T.PPID         = DAT.BEHAVIOR.behavioral_data_table.participant_id;
            T.session      = s .* ones(height(T),1);
            T.session_name = repmat(string(names_sess{s}), height(T), 1);
            results_tables_roi{s} = T;
        end

        results_table_roi = vertcat(results_tables_roi{:});
        fprintf('roi table: %d rows x %d columns\n', ...
                height(results_table_roi), width(results_table_roi));
    end

end


%% ========================================================================
% 3. SAVE
% =========================================================================

if dosave

    fprintf('\n\n'); printhdr('WRITING TABLES'); fprintf('\n\n');

    if do_parcel && ~isempty(results_table_parcel)
        f = fullfile(resultsdir, ['extracted_parcels_', parcel_type, '_', ...
                                  scaling_string, '_', results_suffix_parcel, '.csv']);
        writetable(results_table_parcel, f);
        fprintf('wrote %s\n', f);
    end

    if do_roi && ~isempty(results_table_roi)
        f = fullfile(resultsdir, ['extracted_rois_', scaling_string, '_', ...
                                  results_suffix_roi, '.csv']);
        writetable(results_table_roi, f);
        fprintf('wrote %s\n', f);
    end

end

fprintf('\n\n'); printhdr('DONE'); fprintf('\n\n');
