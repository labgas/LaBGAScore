%% LaBGAScore_mrs_run_osprey_Philips.m
%
%
% *USAGE*
%
% This script runs Osprey analysis of mrs data organized according to
% standard LaBGAS file organization by calling an Osprey jobfile
% (LaBGAScore_mrs_s0_osprey_jobfile) - one jobfile per voxel.
%
% As is the standard, it should be run for the rootdir of the superdataset.
%
% The derivatives & secondlevel subdatasets should be created using datalad commands
% prior to running this script. If this is not the case, it will error out.
%
% Preprocessing of structural (and functional) MRI scans needs to be done
% using fmriprep before running this script, as it uses the segmentations
%
%
% *NOTES*
% 
% 1. example of the BIDS structure created by LaBGAScore_prep_mrs2bids:
%   > BIDS
%       > sub-01
%           > ses-01              
%               > anat
%               > mrs
%                   > sub-01_ses-01_acq-[vox1]press_svs.sdat
%                   > sub-01_ses-01_acq-[vox1]press_svs.spar
%                   > sub-01_ses-01_acq-[vox1]press_ref.sdat
%                   > sub-01_ses-01_acq-[vox1]press_ref.spar
%                   [same 4 files for each voxel]
% 
% 2. based on MRS BIDS proposal as is on 19/02/2023, retrieved via https://docs.google.com/document/d/1pWCb02YNv5W-UZZja24fZrdXLm4X7knXMiZI7E2z7mY/edit#heading=h.ojx3uhh9wmeu
%   (Please keep in mind that final BIDS would require .nii.gz + .json files
%   for MRS data).
%
% 3. If you have a dataset where part of the subjects were scanned on
%       Philips and another part on GE (for example MR8 and PET/MR in
%       LaBGAS case), two different scripts need to be used to
%       analyze those subgroups, calling job files with suffix _Philips and 
%       _GE, respectively, since the eddy current correction option opts.ECC.raw 
%       needs to be set differently for both scanners
%
%
% *OPTION*
%
% whole_analysis = true/false               true runs entire analysis at once, false runs step by step
%
% scanner = 'Philips'/'GE'                  brand of scanner on which data were acquired
%
%
% -------------------------------------------------------------------------
%
% based on script by Melina Hehl
%
% adapted to LaBGAS file organization and looped over voxels by Lukas Van Oudenhove
%
% included option to define scanner type and grab the correct jobfile
%
% date: KU Leuven, February, 2024
%
% -------------------------------------------------------------------------
%
% LaBGAScore_mrs_run_osprey_Philips.m                        v1.4
%
% last modified: 2024/11/20
%
%
%% SET OPTIONS & VARIABLES
% -------------------------------------------------------------------------

whole_analysis = true;

scanner = 'Philips';

study_prefix = 'cfs';                           % prefix used for all scripts in this study

voxelnames = {'pACC'};                          % cell array with voxel names IN THE SAME ORDER AS THEY WERE ACQUIRED!

acq_type = 'press';                             % type of MRS sequence (will be used in 'acq-' label as 'acq-[vox1][Acq_type]')

quant_method = 'TissCorrWaterScaled';           % quantification method

nr_sess = 1;                                    % number of sessions in experiment


%% GET PATHS
% -------------------------------------------------------------------------

% STANDARD PATHS

eval([study_prefix '_prep_s0_define_directories']);

derivrootdir = fullfile(rootdir,'derivatives');
    if ~exist(derivrootdir,'dir')
        error('\nderivatives subdataset "%s" does not exist, please create it using the "datalad create" command and preprocess your MRI data prior to proceeding\n', derivrootdir);
    end

% SANITY CHECK    
    
sourcelist = dir(fullfile(sourcedir,'sub-*'));
sourcesubjs = cellstr(char(sourcelist(:).name));
BIDSlist = dir(fullfile(BIDSdir,'sub-*'));
BIDSsubjs = cellstr(char(BIDSlist(:).name));
derivlist = dir(fullfile(derivdir,'sub-*'));
derivlist = derivlist([derivlist(:).isdir]);
derivsubjs = cellstr(char(derivlist.name));

if isequal(sourcesubjs,BIDSsubjs,derivsubjs)
    warning('\nnumbers and names of subjects in %s, %s, and %s match - good to go',sourcedir,BIDSdir,derivdir);
else
    warning('\nnumbers and names of subjects in %s, %s, and %s do not match - please check before proceeding and make sure all your subjects in the BIDS subdataset are preprocessed',sourcedir,BIDSdir,derivdir);
end


% MRS-SPECIFIC PATHS & FILES

sessions = ones(size(BIDSlist,1),1);

if nr_sess > 1

    for s = 2:nr_sess
        sessions = [sessions; s.*ones(size(BIDSlist,1),1)];
    end
    
end

table_stat = table();
table_stat.subject = repmat(char(BIDSsubjs{:}),nr_sess,1);
table_stat.session = sessions;
table_stat = sortrows(table_stat,'subject');

idx = ones(height(table_stat),1);

for i = 1:size(BIDSlist,1)
    for j = 1:nr_sess
        if nr_sess > 1
            mrsfiles = dir(fullfile(BIDSlist(i).folder,BIDSlist(i).name,['ses-0' num2str(j)],'mrs'));
        else
            mrsfiles = dir(fullfile(BIDSlist(i).folder,BIDSlist(i).name,'mrs')); 
        end
        mrsfiles(1:2)  = [];
        if isempty(mrsfiles) || contains(mrsfiles(1).name,'.7')
            idx(((i-1)*nr_sess)+j) = 0;
        end
        clear mrsfiles
    end
end

table_stat = table_stat(logical(idx),:);


derivospreydir = fullfile(derivrootdir,'osprey');
    if ~exist(derivospreydir,'dir')
        mkdir(derivospreydir);
    end
    
    derivospreyvoxeldirs = cell(size(voxelnames,1));
    files_stat = cell(size(voxelnames,1));
    
    for d = 1:size(voxelnames,1)
        
        derivospreyvoxeldirs{d} = fullfile(derivospreydir,voxelnames{d});
            if ~exist(derivospreyvoxeldirs{d},'dir')
                mkdir(derivospreyvoxeldirs{d});
            end
            
        files_stat{d} = fullfile(derivospreyvoxeldirs{d}, ['stat_' scanner '.csv']);
        writetable(table_stat,files_stat{d});
            
    end
        
    
secondleveldir = fullfile(rootdir,'secondlevel');
    if ~exist(secondleveldir,'dir')
        error('\nsecondlevel subdataset "%s" does not exist, please create it using the "datalad create" command prior to proceeding\n', derivrootdir);
    end

secondlevelmrsdir = fullfile(secondleveldir,'mrs');
    if ~exist(secondlevelmrsdir,'dir')
        mkdir(secondlevelmrsdir);
    end
    
cd(fullfile(derivrootdir,'fmriprep'));

! git annex unannex sub-*/anat/*res-2_label-*_probseg.nii.gz

cd(BIDSdir);

if nr_sess > 1

    ! git annex unannex sub-*/ses-*/anat/*_T1w.nii.gz

else 
    
    ! git annex unannex sub-*/anat/*_T1w.nii.gz
    
end

cd(rootdir);


%% RUN ANALYSIS & EXPORT RESULTS TO EXCEL FILE
% -------------------------------------------------------------------------

% RUN ANALYSIS

for v = 1:size(voxelnames,1)
            
        if nr_sess > 1

            jobFileLocation = fullfile(codedir,'mrs',[study_prefix '_mrs_s0_osprey_jobfile_' voxelnames{v} '_' scanner '.m']);

        else

            jobFileLocation = fullfile(codedir,'mrs',[study_prefix '_mrs_s0_osprey_single_sess_jobfile_' voxelnames{v} '_' scanner '.m']);

        end
        
        
        if whole_analysis

            MRSCont = RunOspreyJob(jobFileLocation); 

        else

            MRSCont = OspreyJob(jobFileLocation); % Generate MRS container
            MRSCont = OspreyLoad(MRSCont); % Load data into MRS container
            MRSCont = OspreyProcess(MRSCont); % Process data
            MRSCont = OspreyFit(MRSCont); % Linear-combination modelling
            MRSCont = OspreyCoreg(MRSCont); % Coregistration module
            MRSCont = OspreySeg(MRSCont); % Segmentation module
            MRSCont = OspreyQuantify(MRSCont); % Quantification module
            MRSCont = OspreyOverview(MRSCont); % Create visualization

        end

    clear jobFileLocation MRSCont
        

% PATHS TO OSPREY RESULTS TSV FILES

dir_OspreyResults = fullfile(derivospreydir, voxelnames{v}, scanner);
path_subjects   = fullfile(dir_OspreyResults, 'subject_names_and_excluded.tsv');
path_quantify   = fullfile(dir_OspreyResults, 'QuantifyResults', 'A_tCr_Voxel_1_Basis_1.tsv');
path_QM         = fullfile(dir_OspreyResults, 'QM_processed_spectra.tsv');


% PATH TO XLS FILE

secondlevelmrsvoxeldir = fullfile(secondlevelmrsdir,voxelnames{v});
    if ~exist(secondlevelmrsvoxeldir,'dir')
        mkdir(secondlevelmrsvoxeldir);
    end

path_result = fullfile(secondlevelmrsvoxeldir, ['OspreyTestResults_' voxelnames{v} '_' scanner '.xlsx']);


% READ RESULTS TSV FILES

results_subjects    = readtable(path_subjects, "FileType","delimitedtext");
results_quantify    = readtable(path_quantify, "FileType","delimitedtext");
results_QM          = readtable(path_QM, "FileType","delimitedtext");
results_QM.Properties.VariableNames = "QC." + results_QM.Properties.VariableNames; % MH add term "QC" to the results of the quality control


% JOIN RESULTS TSV FILES AND WRITE TO EXCEL

results_all = [results_subjects results_QM results_quantify];
results_all.QuantificationMethod = repmat(quant_method, size(results_all, 1), 1); % MH 2024-02-01
writetable(results_all, path_result);

end

