%% LaBGAScore_pet_prep_1_dat_behavioral_data.m
% 
%
% USAGE
%
% This script sets up the DAT structure for second level analysis of PET
% data using CANlab tools in a single PET session design
%
%
% CANLAB NOTES
%
% Optional: Run these load and attach behavioral data from files (e.g., from Excel)            
%
% This script is an example script only. You should modify it to fit your
% needs, which will depend on which types of behavioral/non-imaging data
% you have and what variables you want to store and analyze. 
% 
% The basics are desribed here:
%
% - Store behavioral data tables in any ad hoc format in DAT.BEHAVIOR.
%   This can be a useful reference if you want to change/add custom analyses
%   relating brain to behavior. You can create custom scripts that pull data 
%   from .BEHAVIOR and use it in analyses.
%
% - Store a between-person grouping variable (e.g., patient vs. control,
%   etc.) in DAT.BETWEENPERSON.group. This should be coded with values of 1
%   and -1. Also add fields (see below) for names and colors associated with
%   each group, and a description of what the 1/-1 codes mean.
%   Some analyses consider this variable and run between-group contrasts
%   and/or control for them in analyses of the entire sample (e.g.,
%   "signature response" analyses).  
%   SVM analyses can also be run that use the .group variable. See:
%   prep_3d_run_SVM_betweenperson_contrasts and 
%   c2b_SVM_betweenperson_contrasts  
%
% - If you have no binary group variable,  it is OK to leave the .group
%   field empty. 
% 
% - You can run this script as part of a workflow (prep_1...
%   prep_2...prep_3 etc)
%   You can also run the script AFTER you've prepped all the imaging data, just
%   to add behavioral data to the existing DAT structure.  If so, make sure
%   you RELOAD the existing DAT structure with b_reload_saved_matfiles.m
%   before you run this script.  Otherwise, if you create a new DAT
%   structure, important information saved during the data prep (prep_2...,
%   prep_3...) process will be missing, and you will need to re-run the whole
%   prep sequence.
%
% LABGAS NOTES
%
% - Always make a study-specific copy of this script in your code subdataset, do NOT edit in the repo!
% - This script is highly study-specific, but since we use a fixed
%       BIDS-compatible directory structure and format for the phenotype
%       files, there should be common elements across studies
% - This is an example from a design with a single PET session 
%
%__________________________________________________________________________
%
% modified for pet data by: Lukas Van Oudenhove
% date:   KU Leuven, June 2024
%
%__________________________________________________________________________
% @(#)% LaBGAScore_pet_prep_1_dat_behavioral_data.m         v1.0
% last modified: 2024/06/05
%
%
%% GET PATHS AND OPTIONS
% -------------------------------------------------------------------------

LaBGAScore_pet_a_set_up_paths_always_run_first;

% NOTE: CHANGE THIS TO THE MODEL-SPECIFIC VERSION OF THIS SCRIPT


%% INITIATE DAT STRUCTURE
% -------------------------------------------------------------------------

DAT = struct();
DAT.conditions{1} = modelname;
DAT.functional_wildcard = {'*/*DV_Logan.nii'};    % string to identify firstlevel images


%% READ BEHAVIORAL DATA FROM TSV FILES IN BIDS/PHENOTYPE DIR
% -------------------------------------------------------------------------

phenodir = fullfile(BIDSdir,'phenotype');

% LOAD BEHAVIORAL DATA

behavioral_data_filename = 'ratings_means.tsv';
behavioral_fname_path = fullfile(phenodir, behavioral_data_filename);

if ~exist(behavioral_fname_path, 'file') 
    fprintf(1, 'CANNOT FIND FILE: %s\n',behavioral_fname_path); 
end

behavioral_data_table = readtable(behavioral_fname_path,'TreatAsEmpty','n/a','FileType','text','Delimiter','tab'); % read .tsv file into Matlab table variable

% ADD RATING TABLE TO DAT

DAT.BEHAVIOR = [];

DAT.BEHAVIOR.behavioral_data_table = behavioral_data_table;


%% INITIALIZE GROUP VARIABLE
% -------------------------------------------------------------------------

% INITIALIZE BETWEENPERSON FIELD

DAT.BETWEENPERSON = [];


% INITIALIZE GROUP VARIABLE

% These fields are mandatory, but they can be empty

DAT.BETWEENPERSON.group = [];
DAT.BETWEENPERSON.groupnames = {};
DAT.BETWEENPERSON.groupcolors = {};


%% CHECK DAT, PRINT WARNINGS, SAVE DAT STRUCTURE
% -------------------------------------------------------------------------

printhdr('Save DAT structure and directory names in image_names_and_setup.mat');

savefilename = fullfile(resultsdir, 'image_names_and_setup.mat');
save(savefilename, 'dashes','printstr','printhdr','DAT', 'basedir', 'datadir', 'maskdir', 'resultsdir', 'scriptsdir', 'figsavedir', 'htmlsavedir', '-v7.3');

