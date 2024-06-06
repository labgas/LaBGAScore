%% LaBGAScore_pet_prep_4_apply_signatures_and_save
%
%
% USAGE
%
% This script 
% 1) calculates selected signature responses for PET images
%       included in DAT, and saves them to new fields in DAT.
% 2) calculates responses per NPS subregion if NPS is selected, can/should
%       later be expanded to all signatures?
%
%
% OUTPUT
%
% These fields contain data tables whose columns are PET conditions or contrasts, 
% with variable names based on DAT.conditions or DAT.contrastnames, 
% but with spaces replaced with underscores:
% DAT.SIG_conditions.(myscaling_sigs).(similarity_metric_sigs).(keyword_sigs).signature_name
% DAT.SIG_contrasts.(myscaling_sigs).(similarity_metric_sigs).(keyword_sigs).signature_name
% 
%
% OPTIONS
%
% NOTE: 
% defaults are specified in LaBGAScore_pet_a2_set_default_options, but can be changed below
% in case you want to add signature responses calculated with different
% options, e.g. scaling, and save your new version of the script with a
% letter index
%
% myscaling_sigs = 'raw'/'scaled';
% similarity_metric_sigs = 'dotproduct/'cosine_similarity','correlation';
% keyword_sigs = 'all'/any option from load_image_set;
% 
%
%__________________________________________________________________________
%
% author: Lukas Van Oudenhove
% date:   Leuven, June, 2024
%
%__________________________________________________________________________
% @(#)% LaBGAScore_pet_prep_4_apply_signatures_and_save.m         v1.0
% last modified: 2024/06/06
%
%
%% GET PATHS AND OPTIONS
% -------------------------------------------------------------------------

LaBGAScore_pet_a_set_up_paths_always_run_first;

% NOTE: CHANGE THIS TO THE MODEL-SPECIFIC VERSION OF THIS SCRIPT


%% APPLY SELECTED SIGNATURES
% -------------------------------------------------------------------------

for sig = 1:size(keyword_sigs,2)

    fprintf('\n\n');

    [~,signame] = fileparts(char(keyword_sigs{sig}));
    fprintf(['APPLYING SIGNATURE ', upper(signame), ' ON ', upper(modelname), ' DISTRIBUTION VOLUME ', similarity_metric_sigs]);

    fprintf('\n\n');
    
    if contains(keyword_sigs{sig},filesep) % path to image rather than keyword

        [DAT.SIG.(similarity_metric_sigs).(signame),~] = apply_all_signatures(DATA_OBJ{1}, 'conditionnames', {[modelname '_DV']}, 'similarity_metric', similarity_metric_sigs, 'image_set', keyword_sigs(sig));
        [DAT.SIG.(similarity_metric_sigs).(signame),~] = apply_all_signatures(DATA_OBJ{1}, 'conditionnames', {[modelname '_DV']}, 'similarity_metric', similarity_metric_sigs, 'image_set', keyword_sigs(sig));
        
    else
        
        [DAT.SIG.(similarity_metric_sigs).(signame),~] = apply_all_signatures(DATA_OBJ{1}, 'conditionnames', {[modelname '_DV']}, 'similarity_metric', similarity_metric_sigs, 'image_set', keyword_sigs{sig});
        [DAT.SIG.(similarity_metric_sigs).(signame),~] = apply_all_signatures(DATA_OBJ{1}, 'conditionnames', {[modelname '_DV']}, 'similarity_metric', similarity_metric_sigs, 'image_set', keyword_sigs{sig});
        
    end

end


%% NPS SUBREGIONS
% -------------------------------------------------------------------------

if sum(contains(keyword_sigs,'nps')) > 0 || isequal(keyword_sigs{1},'all')

    % SUBREGION NAMES
    
    posnames = {'vermis' 'rIns' 'rV1' 'rThal' 'lIns' 'rdpIns' 'rS2_Op' 'dACC'};
    negnames = {'rLOC' 'lLOC' 'rpLOC' 'pgACC' 'lSTS' 'rIPL' 'PCC'};

    DAT.SIG.NPSsubregions.posnames = posnames;
    DAT.SIG.NPSsubregions.negnames = negnames;
    
    fprintf('\n\n');
    fprintf('Extracting NPS Subregions, adding to SIG.NPSsubregions');
    fprintf('\n\n');

    % EXTRACT NPS SUBREGIONS

    switch similarity_metric_sigs

        case 'dotproduct'

            [~, ~, ~, DAT.SIG.NPSsubregions.npspos_by_region, DAT.SIG.NPSsubregions.npsneg_by_region] = apply_nps(DATA_OBJ{1}, 'noverbose', 'notables');

        case 'cosine_similatiry'

            [~, ~, ~, DAT.SIG.NPSsubregions.npspos_by_region, DAT.SIG.NPSsubregions.npsneg_by_region] = apply_nps(DATA_OBJ{1}, 'noverbose', 'notables', similarity_metric_sigs);

    end


    % GET AVERAGES FOR NPS SUBREGIONS
   
    DAT.SIG.NPSsubregions.posdat = nanmean([DAT.SIG.NPSsubregions.npspos_by_region{:}])'; % mean across subjects
    DAT.SIG.NPSsubregions.stepos = ste([DAT.SIG.NPSsubregions.npspos_by_region{:}])'; % ste

    DAT.SIG.NPSsubregions.negdat = nanmean([DAT.SIG.NPSsubregions.npsneg_by_region{:}])'; % mean across subjects
    DAT.SIG.NPSsubregions.steneg = ste([DAT.SIG.NPSsubregions.npsneg_by_region{:}])'; % ste

end


%% SAVE RESULTS
% -------------------------------------------------------------------------

fprintf('\n\n');
printhdr('SAVE UPDATED DAT STRUCTURE IN images_names_and_setup.mat');
fprintf('\n\n');

cd(resultsdir); % unannex image_names_and_setup.mat file if already datalad saved to prevent write permission problems
! git annex unannex image_names_and_setup.mat
cd(rootdir);

savefilename = fullfile(resultsdir, 'image_names_and_setup.mat');
save(savefilename, '-append', 'DAT');

