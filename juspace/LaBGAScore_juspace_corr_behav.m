%% LaBGAScore_juspace_corr_behav.m
%
%
% *USAGE*
%
% This script calculates correlations between behavioral responses to conditions and/or contrasts 
% and the results from spatial correlation analysis with PET receptor maps obtained with the JuSpace toolbox
%
% It uses the results file from JuSpace analysis options #4 and #6 for
% conditions and contrasts, respectively
%
% The names of the results files need to start with, respectively
%
% # Results_corrtoPET_4_each_image
%
% # Results_corrtoPET_6_each_image
%
%
% For more info about JuSpace, see the following resources
%
% # JuSpace Github repo: https://github.com/juryxy/JuSpace
%
% # JuSpace HBM paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7814756/
%
%
% *OPTIONS*
%
% # mygroupfieldname: 'contrasts' or 'conditions'
%
% # covs2use: variable name(s) in DAT.BETWEENPERSON.(mygroupnamefield){:} to be correlated with JuSpace results
%
%
% -------------------------------------------------------------------------
%
% author: Lukas Van Oudenhove
%
% date:   KU Leuven, September, 2023
%
% -------------------------------------------------------------------------
%
% LaBGAScore_juspace_corr_behav.m         v1.0
%
% last modified: 2023/09/14
%
%
%% OPTIONS
% -------------------------------------------------------------------------

mygroupfieldname = 'contrasts';
covs2use = 'delta_rating';


%% PREP WORK
% -------------------------------------------------------------------------

% GET MODEL-SPECIFIC PATHS AND OPTIONS

ery_4a_secondlevel_m6m_s0_a_set_up_paths_always_run_first;


% LOAD DSGN, DAT, AND DATA OBJECTS

if ~exist('DSGN','var') || ~exist('DAT','var')
    
    load(fullfile(resultsdir,'image_names_and_setup.mat'));
    
end

% SET JUSPACE DIRS

juspacedir = fullfile(resultsdir,'juspace');
juspaceresultsdir = fullfile(juspacedir, 'results');


%% CALCULATE CORRELATIONS
% -------------------------------------------------------------------------

switch mygroupfieldname
    
    case 'conditions'
        
        corr_results_conditions = cell(1,size(DAT.conditions,2));
        
        for cond = 1:size(DAT.conditions,2)
            
            conddir = fullfile(juspaceresultsdir,DAT.conditions{cond});
            condmatfile = dir(fullfile(conddir,'Results_corrtoPET_4_each_image*.mat'));
            condresults = load(fullfile(condmatfile.folder,condmatfile.name));
            
            petmapnames = cell(1,size(condresults.filesPET,1)+1);
            petmapnames{1} = '';
            rhos = cell(1,size(condresults.filesPET,1)+1);
            rhos{1} = 'rho';
            pvals = cell(1,size(condresults.filesPET,1)+1);
            pvals{1} = 'p-value';
            
            for p = 1:size(condresults.filesPET,1)
                
                petmap = condresults.filesPET{p};
                [~,petmap] = fileparts(petmap);
                petmapnames{1+p} = petmap;
                
                [rhos{1+p},pvals{1+p}] = corr(condresults.res(:,p),DAT.BETWEENPERSON.conditions{cond}.(covs2use));
                
            end
            
            corr_results_conditions{cond} = [petmapnames;rhos;pvals];
            
        end
        
    case 'contrasts'
        
        corr_results_contrasts = cell(1,size(DAT.contrastnames,2));
        
        for cont = 1:size(DAT.contrastnames,2)
            
            contdir = fullfile(juspaceresultsdir,DAT.contrastnames{cont});
            contmatfile = dir(fullfile(contdir,'Results_corrtoPET_6_each_image*.mat'));
            contresults = load(fullfile(contmatfile.folder,contmatfile.name));
            
            petmapnames = cell(1,size(contresults.filesPET,1)+1);
            petmapnames{1} = '';
            rhos = cell(1,size(contresults.filesPET,1)+1);
            rhos{1} = 'rho';
            pvals = cell(1,size(contresults.filesPET,1)+1);
            pvals{1} = 'p-value';
            
            for p = 1:size(contresults.filesPET,1)
                
                petmap = contresults.filesPET{p};
                [~,petmap] = fileparts(petmap);
                petmapnames{1+p} = petmap;
                
                [rhos{1+p},pvals{1+p}] = corr(contresults.stats.res_ind(:,p),DAT.BETWEENPERSON.contrasts{cont}.(covs2use));
            end
            
            corr_results_contrasts{cont} = [petmapnames;rhos;pvals];
            
        end
        
end
