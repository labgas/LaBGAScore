%% LaBGAScore_secondlevel_MS_mat_pipeline.m
%
%
% *USAGE*
%
% MACS Toolbox: Model Space Pipeline
%
% This script assists in setting up a batch for defining a model space that
% can be viewed and executed in the SPM batch editor. It is particularly
% advantageous when the number of subjects or the number of models in your
% analyses is very large.
%
% This script assumes that you have organized your data in the form of a
% subject-model hierarchy looking like this:
%
%     [stat_dir]\
%     [work_dir]\
%         sub01\
%             mod01\
%             mod02\
%             mod03\
%             ...
%             mod08\
%             mod09\
%             mod10\
%         sub02\
%             mod01\
%             ...
%             mod10\
%         sub03\
%         ...
%         sub23\
%         sub24\
%         sub25\
%
% If this is the case, you can simply enter
% - the statistics directory into "stat_dir",
% - the working directory into "work_dir",
% - the subject folder names into "subj_ids" and
% - the model folder names into "mod_nams" below.
%
%
% *OPTIONS*
%
% * ms_name     model space name, used together with stat_dir to determine where the model space directory/file is written
% * ms_suff     model space suffix, used together with ms_name to distinguish different model spaces/analyses from each other
%
%
% *DEPENDENCIES*
%
% MACS Toolbox (SPM batch editor module `spm.tools.MACS.MA_model_space`)
%
%
% *NOTES*
%
% lukasvo76: we have a model-subject hierarchy in our firstlevel subdataset, i.e.
%
%       firstleveldir/
%           mod01/
%               sub01/
%               sub02/
%               ...
%               subxy/
%           mod02/
%               ...
%           ...
%
% -------------------------------------------------------------------------
%
% modified by: Joram Soch (joram.soch@bccn-berlin.de), BCCN Berlin
%              Lukas Van Oudenhove (lukas.vanoudenhove@kuleuven.be), KU Leuven
%
% date:   Joram Soch first edit 18/08/2017, last edit 11/06/2018 (V1.2/V18)
%         Lukas Van Oudenhove adaptation first edit 28/06/2022
%
% -------------------------------------------------------------------------


%%% Step 0: Study parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load LaBGAS standard directory structure
test_prep_s0_define_directories
firstleveldir = fullfile(rootdir,'firstlevel');

% output directory
stat_dir = fullfile(rootdir,'secondlevel','model_selection');

% list of models
mod_nams = dir(fullfile(firstleveldir,'model_*'));
mod_nams = {mod_nams(:).name};

% list of subjects
subj_ids = derivsubjs';
subj_ids2 = dir(fullfile(firstleveldir,mod_nams{1},'sub-*'));
subj_ids2 = {subj_ids2(:).name};

if ~isequal(subj_ids,subj_ids2)
    error('\nsubjects in %s and %s are not identical, please check before proceeding\n', derivdir,fullfile(firstleveldir,mod_nams{1}));
end

% model space details
ms_name  =  'MS01';
ms_suff  =  'test';

% study dimensions
N = numel(subj_ids);
M = numel(mod_nams);


%%% Step 1: Create model space job %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% working directory
clear job
job.dir{1} = strcat(stat_dir,'/',ms_name,'_',ms_suff,'/');

% assemble SPM.mats
for i = 1:N
    for j = 1:M
        job.models{i}{j}{1} = strcat(firstleveldir,'/',mod_nams{j},'/',subj_ids{i},'/','SPM.mat');
    end;
end;

% assemble GLM names
for j = 1:M
    job.names{j} = mod_nams{j};
end;


%%% Step 2: Execute model space job %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% save batch
clear matlabbatch
mkdir(job.dir{1});
filename = strcat(job.dir{1},'batch.mat');
matlabbatch{1}.spm.tools.MACS.MA_model_space = job;
save(filename,'matlabbatch');

% execute job
MA_model_space(job);

% display message
fprintf('\n');
fprintf('\n-> Thank you! The following files have been created:\n');
fprintf('   - SPM batch: %s.\n', strcat(job.dir{1},'batch.mat'));
fprintf('   - MS.mat file: %s.\n', strcat(job.dir{1},'MS.mat'));
fprintf('\n');