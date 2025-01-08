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
% adapted to work with both Philips and GE file formats by LVO
%
% debugged and added deidentification by LVO
%
% date: KU Leuven, February, 2024
%
% -------------------------------------------------------------------------
%
% LaBGAScore_prep_mrs2bids.m                        v1.3
%
% last modified: 2025/01/08
%
%
%% GET PATHS AND DEFINE VOXEL NAMES
% -------------------------------------------------------------------------

LaBGAScore_prep_s0_define_directories; % STUDY-SPECIFIC

voxelnames = {'pACC'};   % cell array with voxel names IN THE SAME ORDER AS THEY WERE ACQUIRED!

acq_type = 'press';             % type of MRS sequence (will be used in 'acq-' label as 'acq-[vox1][Acq_type]')

nr_sess = 1;                    % number of sessions in experiment


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
        plist = dir(fullfile(mrs_subjsourcedir,'P*.7'));
        
            if size(actlist,1) > size(voxelnames,1)*4 || size(reflist,1) > size(voxelnames,1)*4
                error('\nnumber of MRS files should be the same as number of voxelnames (times 4) for subject #%d\n', sub)
            end
            
            if size(plist,1) > size(voxelnames,1)
                error('\nnumber of MRS files should be the same as number of voxelnames for subject #%d\n', sub)
            end
            
            if size(actlist,1) == 0 && size(plist,1) == 0                
               warning('\nno MRS files found for subject #%d, skipping this subject\n', sub) 
               continue
            end
        
            for m = 1:size(voxelnames,1)
                
                voxelname = voxelnames{m};
                
                if size(actlist,1) > 0 % Philips data
                    
                    voxelcounter = m-1;
                
                    actlist_voxel = actlist(1+voxelcounter*2:2+voxelcounter*2,:);
                    reflist_voxel = reflist(1+voxelcounter*2:2+voxelcounter*2,:);
                    
                    cd(mrs_subjsourcedir);
                    PhilipsDeIdentify;
                    cd(rootdir);
                
                    for n = 1:size(actlist_voxel,1)
                        sourceactname = char(actlist_voxel(n).name);
                        sourceactnameparts = strsplit(sourceactname,'.');
                        sourceactext = sourceactnameparts{end};
                        sourceactname_noID = [sourceactnameparts{1} '_noID.' sourceactext];
                        BIDSactname = char(strcat(sourcesubjs{sub},'_acq-',voxelname,acq_type,'_svs.',sourceactext));
                        if ~exist(fullfile(mrs_subjBIDSdir,BIDSactname),'file')
                            movefile(fullfile(mrs_subjsourcedir,sourceactname_noID),fullfile(mrs_subjBIDSdir,BIDSactname));
                        end
                    end

                    for o = 1:size(reflist_voxel,1)
                        sourcerefname = char(reflist_voxel(o).name);
                        sourcerefnameparts = strsplit(sourcerefname,'.');
                        sourcerefext = sourcerefnameparts{end};
                        sourcerefname_noID = [sourcerefnameparts{1} '_noID.' sourcerefext];
                        BIDSrefname = char(strcat(sourcesubjs{sub},'_acq-',voxelname,acq_type,'_ref.',sourcerefext));
                        if ~exist(fullfile(mrs_subjBIDSdir,BIDSrefname),'file')
                            movefile(fullfile(mrs_subjsourcedir,sourcerefname_noID),fullfile(mrs_subjBIDSdir,BIDSrefname));
                        end
                    end
                    
                    clear voxelname actlist_voxel reflist_voxel
                    
                elseif size(plist,1) > 0 % GE data
                    
                    voxelcounter = m;
                    plist_voxel = plist(m);
                    
                    cd(mrs_subjsourcedir);
                    GEDeIdentify;
                    cd(rootdir);
                    
                    for p = 1:size(plist,1)
                        sourcepname = char(plist(p).name);
                        sourcepnameparts = strsplit(sourcepname,'.');
                        sourcepext = sourcepnameparts{end};
                        sourcepname_noID = [sourcepnameparts{1} '_noID.' sourcepext];
                        BIDSpname = char(strcat(sourcesubjs{sub},'_acq-',voxelname,acq_type,'.',sourcepext));
                        if ~exist(fullfile(mrs_subjBIDSdir,BIDSpname),'file')
                            movefile(fullfile(mrs_subjsourcedir,sourcepname_noID),fullfile(mrs_subjBIDSdir,BIDSpname));
                        end
                    end
                    
                    clear voxelname
                   
                end         
                
            end % for loop voxels
            
    elseif nr_sess > 1
        
        for sess = 1:nr_sess
            
            sessid = sprintf('ses-0%d',sess);
            
            mrs_subjsourcedir = fullfile(sourcesubjdirs{sub},sessid,'mrs');
            mrs_subjBIDSdir = fullfile(BIDSsubjdirs{sub},sessid,'mrs');
            actlist = dir(fullfile(mrs_subjsourcedir,'*_act.*'));
            reflist = dir(fullfile(mrs_subjsourcedir,'*_ref.*'));
            plist = dir(fullfile(mrs_subjsourecedir','P*.7'));
        
            if size(actlist,1) > size(voxelnames,1)*4 || size(reflist,1) > size(voxelnames,1)*4
                error('\nnumber of MRS files should be the same as number of voxelnames (times 4) for subject #%d session #%d\n', sub, sess)
            end
            
            if size(plist,2) > size(voxelnames,1)
                error('\nnumber of MRS files should be the same as number of voxelnames for subject #%d session #d\n', sub, sess)
            end
            
            if size(actlist,2) == 0 && size(plist,2) == 0                
               warning('\nno MRS files found for subject #%d session #d, skipping this session\n', sub, sess) 
               continue
            end
        
            for m = 1:size(voxelnames,1)
                
                voxelname = voxelnames{m};
                
                if size(actlist,1) > 0 % Philips data
                
                    voxelcounter = m-1;

                    actlist_voxel = actlist(1+voxelcounter*2:2+voxelcounter*2,:);
                    reflist_voxel = reflist(1+voxelcounter*2:2+voxelcounter*2,:);
                    
                    cd(mrs_subjsourcedir);
                    PhilipsDeIdentify;
                    cd(rootdir);
                
                    for n = 1:size(actlist_voxel,1)
                        sourceactname = char(actlist_voxel(n).name);
                        sourceactnameparts = strsplit(sourceactname,'.');
                        sourceactext = sourceactnameparts{end};
                        sourceactname_noID = [sourceactnameparts{1} '_noID.' sourceactext];
                        BIDSactname = char(strcat(sourcesubjs{sub},'_',sessid,'_acq-',voxelname,acq_type,'_svs.',sourceactext));
                        if ~exist(fullfile(mrs_subjBIDSdir,BIDSactname),'file')
                            movefile(fullfile(mrs_subjsourcedir,sourceactname_noID),fullfile(mrs_subjBIDSdir,BIDSactname));
                        end
                    end

                    for o = 1:size(reflist_voxel,1)
                        sourcerefname = char(reflist_voxel(o).name);
                        sourcerefnameparts = strsplit(sourcerefname,'.');
                        sourcerefext = sourcerefnameparts{end};
                        sourcerefname_noID = [sourcerefnameparts{1} '_noID.' sourcerefext];
                        BIDSrefname = char(strcat(sourcesubjs{sub},'_',sessid,'_acq-',voxelname,acq_type,'_ref.',sourcerefext));
                        if ~exist(fullfile(mrs_subjBIDSdir,BIDSrefname),'file')
                            movefile(fullfile(mrs_subjsourcedir,sourcerefname_noID),fullfile(mrs_subjBIDSdir,BIDSrefname));
                        end
                    end
                    
                    clear voxelname actlist_voxel reflist_voxel
                    
                elseif size(plist,1) > 0 % GE data
                    
                    voxelcounter = m;
                    plist_voxel = plist(m);
                    
                    cd(mrs_subjsourcedir);
                    GEDeIdentify;
                    cd(rootdir);
                    
                    for p = 1:size(plist,1)
                        sourcepname = char(plist(n).name);
                        sourcepnameparts = strsplit(sourcepname,'.');
                        sourcepext = sourcepnameparts{end};
                        sourcepname_noID = [sourcepnameparts{1} '_noID.' sourcepext];
                        BIDSpname = char(strcat(sourcesubjs{sub},'_',sessid,'_acq-',voxelname,acq_type,'.',sourcepext));
                        if ~exist(fullfile(mrs_subjBIDSdir,BIDSpname),'file')
                            movefile(fullfile(mrs_subjsourcedir,sourcepname_noID),fullfile(mrs_subjBIDSdir,BIDSpname));
                        end
                    end
                    
                    clear voxelname
                   
                end
                
            end % for loop voxels
            
            clear sessid mrs_subjsourcedir mrs_subjBIDSdir actlist reflist

        end % for loop sessions
               
    end % if loop sessions
    
end % for loop over subjects

