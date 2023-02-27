%%% LaBGAScore_prep_parrec2bids
%
% This script converts Philips PARREC files to Nifti files and names and
% organizes them according to BIDS specification
%
% USAGE
%
% Script should be run from the root directory of the superdataset, in this
% example case /data/proj_bitter-reward, and assumes that 
% 1. the sourcedata, BIDS, and code subdataset are already created (using datalad if you use it)
% 2. the sourcedata subdataset is organized according to LaBGAS convention
% - see fmri analysis workflow Google doc on LaBGAS drive
%
% The script is study-specific, I indicate in the code below where
% study-specific changes will need to be made
%
% DEPENDENCIES
% 
% 1. spm on your Matlab path
%       https://www.fil.ion.ucl.ac.uk/spm/
% 2. Xiangruili's dicm2nii Github repo on your Matlab path
%       https://github.com/xiangruili/dicm2nii
%
%__________________________________________________________________________
%
% author: lukas.vanoudenhove@kuleuven.be
% date: August, 2022
%
%__________________________________________________________________________
% @(#)% LaBGAScore_prep_parrec2bids.m           v1.2
% last modified: 2022/12/07
%
%
%% STUDY-SPECIFIC OPTIONS AND SETTINGS
%--------------------------------------------------------------------------

nr_sess = 2; % number of sessions for each subject
task1name = 'food_images'; % names of task(s)
task2name = 'FID';
slice_scan_order_exam_card = 'FH'; % get this info from your exam card, other options are 'HF', and 'interleaved'
fold_over_direction_exam_card = 'AP'; % other options unknown, but this is standard on MR8
fat_shift_direction_exam_card = 'P'; % other option 'A';


%% DETERMINE PHASEENCODINGDIRECTION FOR JSON FILES FROM INFO ABOVE
%--------------------------------------------------------------------------

if strcmpi(fold_over_direction_exam_card,'AP')
    if strcmpi(fat_shift_direction_exam_card,'P')
        phase_encoding_string = 'j';
    elseif strcmpi(fat_shift_direction_exam_card,'A')
        phase_encoding_string = 'j-';
    else
        error('\nWrong option "%s" specified in fat_shift_direction_exam_card variable, check your exam card and choose between "A", or "P"\n',fat_shift_direction_exam_card)
    end

else
    error('\nUnknown option "%s" specified in fold_over_direction_exam_card variable, currently only "AP" is implemented, check with your MR physicist before proceeding\n',fold_over_direction_exam_card)

end
                    
                    
%% DEFINE DIRECTORIES AND ADD CODE DIR TO MATLAB PATH
%--------------------------------------------------------------------------

rootdir = pwd;
githubrootdir = '/data/master_github_repos';
sourcedir = fullfile(rootdir,'sourcedata');
BIDSdir = fullfile(rootdir,'BIDS');
codedir = fullfile(rootdir,'code');

matlabpath = path;

if ~exist('spm.m','file')
    spmpathcommand = "addpath('your_spm_rootdir','-end')";
    error('\nspm12 not found on Matlab path, please add WITHOUT subfolders using the Matlab GUI or type %s in Matlab terminal before proceeding',spmpathcommand)

else
    spmrootdir = which('spm.m');
    [spmrootdir,~,~] = fileparts(spmrootdir);

end
    
if sum(contains(matlabpath,codedir)) == 0
    addpath(genpath(codedir),'-end');
    warning('\nadding %s to end of Matlab path\n',codedir)

end


%% READ IN SUBJECT LIST FROM SOURCEDATA AND WRITE IN BIDSDIR
%--------------------------------------------------------------------------

sourcelist = dir(fullfile(sourcedir,'sub-*'));
sourcesubjs = cellstr(char(sourcelist(:).name));
sourcesubjdirs = cell(size(sourcesubjs,1),1);

for sourcesub = 1:size(sourcesubjs,1)    
    
    sourcesubjdirs{sourcesub,1} = fullfile(sourcelist(sourcesub).folder, sourcelist(sourcesub).name);

end

for sub = 1:size(sourcesubjdirs,1)
    
    subjBIDSdir = fullfile(BIDSdir,sourcesubjs{sub});
    
    if ~exist(subjBIDSdir,'dir')
        mkdir(subjBIDSdir);
    end
    
        if nr_sess > 1
            
            for sess = 1:nr_sess
                sessid = sprintf('ses-%d',sess);
                subjsessdir = fullfile(subjBIDSdir,sessid);
                if ~exist(subjsessdir,'dir')
                    mkdir(subjsessdir);
                    mkdir(subjsessdir,'anat');
                    mkdir(subjsessdir,'func');
                end
            end
            
        end
        
    clear subjBIDSdir sess sessid subjsessdir
end

clear sub


BIDSlist = dir(fullfile(BIDSdir,'sub-*'));
BIDSsubjs = cellstr(char(BIDSlist(:).name));
BIDSsubjdirs = cell(size(BIDSsubjs,1),1);

for BIDSsub = 1:size(BIDSsubjs,1)    
    BIDSsubjdirs{BIDSsub,1} = fullfile(BIDSlist(BIDSsub).folder, BIDSlist(BIDSsub).name);

end


%% CONVERT PARRECS TO NIFTI AND WRITE TO THE RIGHT PLACE
%--------------------------------------------------------------------------

for sub = 1:size(sourcesubjdirs,1)
    
    if nr_sess == 1
        
        subjsourcedir = fullfile(sourcedir,sourcesubjs{sub},'PARREC');
        subjBIDSdir = fullfile(BIDSdir, sourcesubjs{sub});
        
        dicm2nii(char(subjsourcedir),char(subjBIDSdir),1);
            
        load(fullfile(subjBIDSdir,'dcmHeaders.mat'));

        parrec_headers = fields(h);

        anatniilist = dir(fullfile(subjBIDSdir,'*_3DTFE_*.nii.gz')); % STUDY-SPECIFIC: THE NAMES OF YOUR STRUCTURAL PARRECS MAY DIFFER
            for i = 1:size(anatniilist,1)
                sourceanatniiname = char(anatniilist(i).name);
                BIDSanatniiname = char(strcat(sourcesubjs{sub},'_T1w.nii.gz'));
                movefile(fullfile(subjBIDSdir,sourceanatniiname),fullfile(subjBIDSdir,'anat',BIDSanatniiname));
            end

        anatjsonlist = dir(fullfile(subjBIDSdir,'*_3DTFE_*.json')); % STUDY-SPECIFIC: THE NAMES OF YOUR STRUCTURAL PARRECS MAY DIFFER
            for i = 1:size(anatjsonlist,1)
                sourceanatjsonname = char(anatjsonlist(i).name);
                anatjson = spm_jsonread(fullfile(subjBIDSdir,sourceanatjsonname));
                    anatjson.SeriesDescription = [];
                spm_jsonwrite(fullfile(subjBIDSdir,sourceanatjsonname),anatjson);
                BIDSanatjsonname = char(strcat(sourcesubjs{sub},'_T1w.json'));
                movefile(fullfile(subjBIDSdir,sourceanatjsonname),fullfile(subjBIDSdir,'anat',BIDSanatjsonname));
            end

        task1niilist1 = dir(fullfile(subjBIDSdir,'*_Exp1_*.nii.gz')); % STUDY-SPECIFIC: THE NAMES OF YOUR FUNCTIONAL PARRECS MAY DIFFER
        task1niilist2 = dir(fullfile(subjBIDSdir,'*_Exp2_*.nii.gz'));
        task1niilist = [task1niilist1;task1niilist2];
            for i = 1:size(task1niilist,1)
                runid = sprintf('_run-%d',i);
                sourcetask1niiname = char(task1niilist(i).name);
                BIDStask1niiname = char(strcat(sourcesubjs{sub},'_task-', task1name, runid, '_bold.nii.gz'));
                movefile(fullfile(subjBIDSdir,sourcetask1niiname),fullfile(subjBIDSdir,'func',BIDStask1niiname));
            end

        task1jsonlist1 = dir(fullfile(subjBIDSdir,'*_Exp1_*.json')); % STUDY-SPECIFIC: THE NAMES OF YOUR FUNCTIONAL PARRECS MAY DIFFER
        task1jsonlist2 = dir(fullfile(subjBIDSdir,'*_Exp2_*.json'));
        task1jsonlist = [task1jsonlist1;task1jsonlist2];
            for i = 1:size(task1jsonlist,1)
                runid = sprintf('_run-%d',i);
                sourcetask1jsonname = char(task1jsonlist(i).name);

                % the following is based on Chris Rorden's code 
                % https://neurostars.org/t/deriving-slice-timing-order-from-philips-par-rec-and-console-information/17688

                n_slices = h.(parrec_headers{size(anatjsonlist,1)+i}).LocationsInAcquisition;
                TR = (h.(parrec_headers{size(anatjsonlist,1)+i}).RepetitionTime)/1000; % convert to seconds
                last_slice_time = TR - TR/n_slices;

                switch slice_scan_order_exam_card
                    case 'FH' % ascending
                        slice_times = linspace(0, last_slice_time, n_slices); 
                    case 'HF' % descending
                        slice_times = flip(linspace(0, last_slice_time, n_slices));
                    case 'interleaved'
                        order = [1:2:n_slices 2:2:n_slices]; 
                        slice_times = linspace(0, last_slice_time, n_slices); 
                        slice_times(order) = slice_times;
                    otherwise
                        error('\nWrong option "%s" specified in slice_scan_order_exam_card variable, check your exam card and choose between "FH", HF", or "interleaved"\n',slice_scan_order_exam_card)
                end

                task1json = spm_jsonread(fullfile(subjBIDSdir,sourcetask1jsonname));
                    task1json.SeriesDescription = [];
                    task1json.PhaseEncodingDirection = phase_encoding_string;
                    task1json.TaskName = task1name;
                    task1json.SliceTiming = slice_times;
                spm_jsonwrite(fullfile(subjBIDSdir,sourcetask1jsonname),task1json);

                BIDStask1jsonname = char(strcat(sourcesubjs{sub},'_task-', task1name , runid, '_bold.json'));

                movefile(fullfile(subjBIDSdir,sourcetask1jsonname),fullfile(subjBIDSdir,'func',BIDStask1jsonname));
            end

        task2niilist1 = dir(fullfile(subjBIDSdir,'*_Exp3_*.nii.gz')); % STUDY-SPECIFIC: THE NAMES OF YOUR FUNCTIONAL PARRECS MAY DIFFER
        task2niilist2 = dir(fullfile(subjBIDSdir,'*_Exp4_*.nii.gz'));
        task2niilist = [task2niilist1;task2niilist2];
            for i = 1:size(task2niilist,1)
                runid = sprintf('_run-%d',i);
                sourcetask2niiname = char(task2niilist(i).name);
                BIDStask2niiname = char(strcat(sourcesubjs{sub},'_task-', task2name, runid, '_bold.nii.gz'));
                movefile(fullfile(subjBIDSdir,sourcetask2niiname),fullfile(subjBIDSdir,'func',BIDStask2niiname));
            end

        task2jsonlist1 = dir(fullfile(subjBIDSdir,'*_Exp3_*.json')); % STUDY-SPECIFIC: THE NAMES OF YOUR FUNCTIONAL PARRECS MAY DIFFER
        task2jsonlist2 = dir(fullfile(subjBIDSdir,'*_Exp4_*.json')); 
        task2jsonlist = [task2jsonlist1;task2jsonlist2];
            for i = 1:size(task2jsonlist,1)
                runid = sprintf('_run-%d',i);
                sourcetask2jsonname = char(task2jsonlist(i).name);

                % the following is based on Chris Rorden's code 
                % https://neurostars.org/t/deriving-slice-timing-order-from-philips-par-rec-and-console-information/17688

                n_slices = h.(parrec_headers{size(anatjsonlist,1)+size(task1jsonlist,1)+i}).LocationsInAcquisition;
                TR = (h.(parrec_headers{size(anatjsonlist,1)+size(task1jsonlist,1)+i}).RepetitionTime)/1000; % convert to seconds
                last_slice_time = TR - TR/n_slices;

                switch slice_scan_order_exam_card
                    case 'FH' % ascending
                        slice_times = linspace(0, last_slice_time, n_slices); 
                    case 'HF' % descending
                        slice_times = flip(linspace(0, last_slice_time, n_slices));
                    case 'interleaved'
                        order = [1:2:n_slices 2:2:n_slices]; 
                        slice_times = linspace(0, last_slice_time, n_slices); 
                        slice_times(order) = slice_times;
                    otherwise
                        error('\nWrong option "%s" specified in slice_scan_order_exam_card variable, check your exam card and choose between "FH", HF", or "interleaved"\n',slice_scan_order_exam_card)
                end

                task2json = spm_jsonread(fullfile(subjBIDSdir,sourcetask2jsonname));
                    task2json.SeriesDescription = [];
                    task2json.PhaseEncodingDirection = phase_encoding_string;
                    task2json.TaskName = task2name;
                    task2json.SliceTiming = slice_times;
                spm_jsonwrite(fullfile(subjBIDSdir,sourcetask2jsonname),task2json);

                BIDStask2jsonname = char(strcat(sourcesubjs{sub},'_task-', task2name, runid, '_bold.json'));

                movefile(fullfile(subjBIDSdir,sourcetask2jsonname),fullfile(subjBIDSdir,'func',BIDStask2jsonname));
            end
        
    else
        
        for sess = 1:nr_sess
            
            sessid = sprintf('ses-%d',sess);
            
            subjsourcesessdir = fullfile(sourcedir,sourcesubjs{sub},sessid,'PARREC');
            subjBIDSsessdir = fullfile(BIDSdir,sourcesubjs{sub},sessid);
            
            dicm2nii(char(subjsourcesessdir),char(subjBIDSsessdir),1);
            
            load(fullfile(subjBIDSsessdir,'dcmHeaders.mat'));

            parrec_headers = fields(h);
            
            anatniilist = dir(fullfile(subjBIDSsessdir,'*_3DTFE_*.nii.gz')); % STUDY-SPECIFIC: THE NAMES OF YOUR STRUCTURAL PARRECS MAY DIFFER
                for i = 1:size(anatniilist,1)
                    sourceanatniiname = char(anatniilist(i).name);
                    BIDSanatniiname = char(strcat(sourcesubjs{sub},'_',sessid,'_T1w.nii.gz'));
                    movefile(fullfile(subjBIDSsessdir,sourceanatniiname),fullfile(subjBIDSsessdir,'anat',BIDSanatniiname));
                end
            
            anatjsonlist = dir(fullfile(subjBIDSsessdir,'*_3DTFE_*.json')); % STUDY-SPECIFIC: THE NAMES OF YOUR STRUCTURAL PARRECS MAY DIFFER
                for i = 1:size(anatjsonlist,1)
                    sourceanatjsonname = char(anatjsonlist(i).name);
                    anatjson = spm_jsonread(fullfile(subjBIDSsessdir,sourceanatjsonname));
                        anatjson.SeriesDescription = [];
                    spm_jsonwrite(fullfile(subjBIDSsessdir,sourceanatjsonname),anatjson);
                    BIDSanatjsonname = char(strcat(sourcesubjs{sub},'_',sessid,'_T1w.json'));
                    movefile(fullfile(subjBIDSsessdir,sourceanatjsonname),fullfile(subjBIDSsessdir,'anat',BIDSanatjsonname));
                end

            task1niilist1 = dir(fullfile(subjBIDSsessdir,'*_Exp1_*.nii.gz')); % STUDY-SPECIFIC: THE NAMES OF YOUR FUNCTIONAL PARRECS MAY DIFFER
            task1niilist2 = dir(fullfile(subjBIDSsessdir,'*_Exp2_*.nii.gz'));
            task1niilist = [task1niilist1;task1niilist2];
                for i = 1:size(task1niilist,1)
                    runid = sprintf('_run-%d',i);
                    sourcetask1niiname = char(task1niilist(i).name);
                    BIDStask1niiname = char(strcat(sourcesubjs{sub},'_',sessid,'_task-', task1name, runid, '_bold.nii.gz'));
                    movefile(fullfile(subjBIDSsessdir,sourcetask1niiname),fullfile(subjBIDSsessdir,'func',BIDStask1niiname));
                end

            task1jsonlist1 = dir(fullfile(subjBIDSsessdir,'*_Exp1_*.json')); % STUDY-SPECIFIC: THE NAMES OF YOUR FUNCTIONAL PARRECS MAY DIFFER
            task1jsonlist2 = dir(fullfile(subjBIDSsessdir,'*_Exp2_*.json'));
            task1jsonlist = [task1jsonlist1;task1jsonlist2];
                for i = 1:size(task1jsonlist,1)
                    runid = sprintf('_run-%d',i);
                    sourcetask1jsonname = char(task1jsonlist(i).name);
                    
                    % the following is based on Chris Rorden's code 
                    % https://neurostars.org/t/deriving-slice-timing-order-from-philips-par-rec-and-console-information/17688
                    
                    n_slices = h.(parrec_headers{size(anatjsonlist,1)+i}).LocationsInAcquisition;
                    TR = (h.(parrec_headers{size(anatjsonlist,1)+i}).RepetitionTime)/1000; % convert to seconds
                    last_slice_time = TR - TR/n_slices;
                    
                    switch slice_scan_order_exam_card
                        case 'FH' % ascending
                            slice_times = linspace(0, last_slice_time, n_slices); 
                        case 'HF' % descending
                            slice_times = flip(linspace(0, last_slice_time, n_slices));
                        case 'interleaved'
                            order = [1:2:n_slices 2:2:n_slices]; 
                            slice_times = linspace(0, last_slice_time, n_slices); 
                            slice_times(order) = slice_times;
                        otherwise
                            error('\nWrong option "%s" specified in slice_scan_order_exam_card variable, check your exam card and choose between "FH", HF", or "interleaved"\n',slice_scan_order_exam_card)
                    end
                    
                    task1json = spm_jsonread(fullfile(subjBIDSsessdir,sourcetask1jsonname));
                        task1json.SeriesDescription = [];
                        task1json.PhaseEncodingDirection = phase_encoding_string;
                        task1json.TaskName = task1name;
                        task1json.SliceTiming = slice_times;
                    spm_jsonwrite(fullfile(subjBIDSsessdir,sourcetask1jsonname),task1json);
                    
                    BIDStask1jsonname = char(strcat(sourcesubjs{sub},'_task-', task1name , runid, '_bold.json'));
                    
                    movefile(fullfile(subjBIDSsessdir,sourcetask1jsonname),fullfile(subjBIDSsessdir,'func',BIDStask1jsonname));
                end

            task2niilist1 = dir(fullfile(subjBIDSsessdir,'*_Exp3_*.nii.gz')); % STUDY-SPECIFIC: THE NAMES OF YOUR FUNCTIONAL PARRECS MAY DIFFER
            task2niilist2 = dir(fullfile(subjBIDSsessdir,'*_Exp4_*.nii.gz'));
            task2niilist = [task2niilist1;task2niilist2];
                for i = 1:size(task2niilist,1)
                    runid = sprintf('_run-%d',i);
                    sourcetask2niiname = char(task2niilist(i).name);
                    BIDStask2niiname = char(strcat(sourcesubjs{sub},'_',sessid,'_task-', task2name, runid, '_bold.nii.gz'));
                    movefile(fullfile(subjBIDSsessdir,sourcetask2niiname),fullfile(subjBIDSsessdir,'func',BIDStask2niiname));
                end

            task2jsonlist1 = dir(fullfile(subjBIDSsessdir,'*_Exp3_*.json')); % STUDY-SPECIFIC: THE NAMES OF YOUR FUNCTIONAL PARRECS MAY DIFFER
            task2jsonlist2 = dir(fullfile(subjBIDSsessdir,'*_Exp4_*.json')); 
            task2jsonlist = [task2jsonlist1;task2jsonlist2];
                for i = 1:size(task2jsonlist,1)
                    runid = sprintf('_run-%d',i);
                    sourcetask2jsonname = char(task2jsonlist(i).name);
                    
                    % the following is based on Chris Rorden's code 
                    % https://neurostars.org/t/deriving-slice-timing-order-from-philips-par-rec-and-console-information/17688
                    
                    n_slices = h.(parrec_headers{size(anatjsonlist,1)+size(task1jsonlist,1)+i}).LocationsInAcquisition;
                    TR = (h.(parrec_headers{size(anatjsonlist,1)+size(task1jsonlist,1)+i}).RepetitionTime)/1000; % convert to seconds
                    last_slice_time = TR - TR/n_slices;
                    
                    switch slice_scan_order_exam_card
                        case 'FH' % ascending
                            slice_times = linspace(0, last_slice_time, n_slices); 
                        case 'HF' % descending
                            slice_times = flip(linspace(0, last_slice_time, n_slices));
                        case 'interleaved'
                            order = [1:2:n_slices 2:2:n_slices]; 
                            slice_times = linspace(0, last_slice_time, n_slices); 
                            slice_times(order) = slice_times;
                        otherwise
                            error('\nWrong option "%s" specified in slice_scan_order_exam_card variable, check your exam card and choose between "FH", HF", or "interleaved"\n',slice_scan_order_exam_card)
                    end
                    
                    task2json = spm_jsonread(fullfile(subjBIDSsessdir,sourcetask2jsonname));
                        task2json.SeriesDescription = [];
                        task2json.PhaseEncodingDirection = phase_encoding_string;
                        task2json.TaskName = task2name;
                        task2json.SliceTiming = slice_times;
                    spm_jsonwrite(fullfile(subjBIDSsessdir,sourcetask2jsonname),task2json);
                    
                    BIDStask2jsonname = char(strcat(sourcesubjs{sub},'_',sessid,'_task-', task2name, runid, '_bold.json'));
                    
                    movefile(fullfile(subjBIDSsessdir,sourcetask2jsonname),fullfile(subjBIDSsessdir,'func',BIDStask2jsonname));
                end
                
                movefile(fullfile(subjBIDSsessdir,'dcmHeaders.mat'), fullfile(subjsourcesessdir,'dcmHeaders.mat')); % move to sourcedata to prevent BIDS validator from erroring
            
        end % for loop over sessions
        
    end % if loop single versus multiple sessions
    
end % for loop over subjects

