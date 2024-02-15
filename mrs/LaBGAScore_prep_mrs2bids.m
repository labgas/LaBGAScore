%% LaBGAScore_prep_mrs2bids.m
%
%
% *USAGE*
%
% This script converts raw mrs data (as exported from Philips scanner and stored 
% in sourcedata folder under standard LaBGAS file organization) to BIDS structure.
%
% As is the standard, it should be run for the rootdir of the superdataset.
%
%
% *NOTES*
% 
% 1. example of the BIDS structure this script creates:
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
%
% -------------------------------------------------------------------------
%
% based on script by Melina Hehl
%
% adapted to LaBGAS file organization by Lukas Van Oudenhove
%
% date: KU Leuven, February, 2024
%
% -------------------------------------------------------------------------
%
% LaBGAScore_prep_mrs2bids.m                        v1.0
%
% last modified: 2024/02/12
%
%
%% GET PATHS AND DEFINE VOXEL NAMES
% -------------------------------------------------------------------------

moodbugs2_prep_s0_define_directories;

voxelnames = {'LINS';'RINS'};   % cell array with voxel names IN THE SAME ORDER AS THEY WERE ACQUIRED!

acq_type = 'press';             % type of MRS sequence (will be used in 'acq-' label as 'acq-[vox1][Acq_type]')

nr_sess = 2;                    % number of sessions in experiment


%% READ IN SUBJECT LIST FROM SOURCEDATA AND WRITE IN BIDSDIR
%--------------------------------------------------------------------------

sourcelist = dir(fullfile(sourcedir,'sub-*'));
sourcesubjs = cellstr(char(sourcelist(:).name));
sourcesubjdirs = cell(size(sourcesubjs,1),1);

if nr_sess == 1
    mrs_sourcesubjdirs = cell(size(sourcesubjs,1),1);
end

for sourcesub = 1:size(sourcesubjs,1)    
    
    sourcesubjdirs{sourcesub,1} = fullfile(sourcelist(sourcesub).folder, sourcelist(sourcesub).name);
    
        if nr_sess == 1
            mrs_sourcesubjdirs{sourcesub,1} = fullfile(sourcelist(sourcesub).folder, sourcelist(sourcesub).name, 'mrs');
        end

end

for sub = 1:size(sourcesubjdirs,1)
    
    subjBIDSdir = fullfile(BIDSdir,sourcesubjs{sub});
    
    if ~exist(subjBIDSdir,'dir')
        mkdir(subjBIDSdir);
    end
    
    if nr_sess == 1
        
        mrs_subjBIDSdir = fullfile(BIDSdir,sourcesubjs{sub},'mrs');
    
        if ~exist(mrs_subjBIDSdir,'dir')
            mkdir(mrs_subjBIDSdir);
        end
    
    elseif nr_sess > 1
            
            for sess = 1:nr_sess
                
                sessid = sprintf('ses-0%d',sess);
                
                subjsessdir = fullfile(subjBIDSdir,sessid);
                mrs_subjsessdir = fullfile(subjBIDSdir,sessid,'mrs');
                
                if ~exist(subjsessdir,'dir')
                    mkdir(subjsessdir);
                end
                
                if ~exist(mrs_subjsessdir,'dir')
                    mkdir(mrs_subjsessdir);
                end
                
            end
            
    end
        
    clear subjBIDSdir mrs_subjBIDSdir sess sessid subjsessdir mrs_subjsessdir
    
end

clear sub


BIDSlist = dir(fullfile(BIDSdir,'sub-*'));
BIDSsubjs = cellstr(char(BIDSlist(:).name));
BIDSsubjdirs = cell(size(BIDSsubjs,1),1);

if nr_sess == 1
    mrs_BIDSsubjdirs = cell(size(BIDSsubjs,1),1);
end

for BIDSsub = 1:size(BIDSsubjs,1)
    
    BIDSsubjdirs{BIDSsub,1} = fullfile(BIDSlist(BIDSsub).folder, BIDSlist(BIDSsub).name);
    
    if nr_sess == 1
        mrs_BIDSsubjdirs{BIDSsub,1} = fullfile(BIDSsubjdirs{BIDSsub,1}, 'mrs');
    end
    
end


%% COPY MRS DATA AND RENAME TO BIDS CONVENTION
% -------------------------------------------------------------------------

for sub = 1:size(sourcesubjdirs,1)
    
    if nr_sess == 1
        
        mrs_subjsourcedir = mrs_sourcesubjdirs{sub};
        mrs_subjBIDSdir = mrs_BIDSsubjdirs{sub};
        actlist = dir(fullfile(mrs_subjsourcedir,'*_act.*'));
        reflist = dir(fullfile(mrs_subjsourcedir,'*_ref.*'));
        
            if size(actlist,2) > size(voxelnames,1)*4 || size(reflist,2) > size(voxelnames,1)*4
                error('\nnumber of MRS files should be the same as number of voxelnames (times 4) for subject #%d\n', sub)
            end
        
            for m = 1:size(voxelnames,1)
                
                voxelname = voxelnames{m};
                
                for n = 1:size(actlist,1)
                    sourceactname = char(actlist(n).name);
                    sourceactnameparts = strsplit(sourceactname,'.');
                    sourceactext = sourceactnameparts{end};
                    BIDSactname = char(strcat(sourcesubjs{sub},'_acq-',voxelname,acq_type,'_svs.',sourceactext));
                    copyfile(fullfile(mrs_subjsourcedir,sourceactname),fullfile(mrs_subjBIDSdir,BIDSactname));
                end
                
                for o = 1:size(reflist,1)
                    sourcerefname = char(reflist(o).name);
                    sourcerefnameparts = strsplit(sourcerefname,'.');
                    sourcerefext = sourcerefnameparts{end};
                    BIDSrefname = char(strcat(sourcesubjs{sub},'_acq-',voxelname,acq_type,'_ref.',sourcerefext));
                    copyfile(fullfile(mrs_subjsourcedir,sourcerefname),fullfile(mrs_subjBIDSdir,BIDSrefname));
                end
                
                clear voxelname
                
            end % for loop voxels
            
    elseif nr_sess > 1
        
        for sess = 1:nr_sess
            
            sessid = sprintf('ses-0%d',sess);
            
            mrs_subjsourcedir = fullfile(sourcesubjdirs{sub},sessid,'mrs');
            mrs_subjBIDSdir = fullfile(BIDSsubjdirs{sub},sessid,'mrs');
            actlist = dir(fullfile(mrs_subjsourcedir,'*_act.*'));
            reflist = dir(fullfile(mrs_subjsourcedir,'*_ref.*'));
        
            if size(actlist,1) > size(voxelnames,1)*4 || size(reflist,1) > size(voxelnames,1)*4
                error('\nnumber of MRS files should be the same as number of voxelnames (time 4) for subject #%d session #%d\n', sub, sess)
            end
        
            for m = 1:size(voxelnames,1)
                
                voxelname = voxelnames{m};
                
                for n = 1:size(actlist,1)
                    sourceactname = char(actlist(n).name);
                    sourceactnameparts = strsplit(sourceactname,'.');
                    sourceactext = sourceactnameparts{end};
                    BIDSactname = char(strcat(sourcesubjs{sub},'_',sessid,'_acq-',voxelname,acq_type,'_svs.',sourceactext));
                    copyfile(fullfile(mrs_subjsourcedir,sourceactname),fullfile(mrs_subjBIDSdir,BIDSactname));
                end
                
                for o = 1:size(reflist,1)
                    sourcerefname = char(reflist(o).name);
                    sourcerefnameparts = strsplit(sourcerefname,'.');
                    sourcerefext = sourcerefnameparts{end};
                    BIDSrefname = char(strcat(sourcesubjs{sub},'_',sessid,'_acq-',voxelname,acq_type,'_ref.',sourcerefext));
                    copyfile(fullfile(mrs_subjsourcedir,sourcerefname),fullfile(mrs_subjBIDSdir,BIDSrefname));
                end
                
                clear voxelname
                
            end % for loop voxels
            
            clear sessid mrs_subjsourcedir mrs_subjBIDSdir actlist reflist

        end % for loop sessions
               
    end % if loop sessions
    
end % for loop over subjects

