%% LaBGAScore_prep_s1_write_events_tsv
%
% This script reads logfiles, extracts the onsets, durations, and ratings for
% different conditions, and writes events.tsv files to the BIDS dir for
% each subject
% It also contains an option to write a single phenotype file with
% trial-by-trial ratings for all subjects
% 
% USAGE
%
% Script should be run from the root directory of the superdataset, e.g.
% /data/proj_discoverie
% The script is highly study-specific, as logfiles will vary with design,
% stimulus presentation software used, etc
% Hence, it is provided in LaBGAScore as an example and needs to be
% downloaded and adapted to the code subdataset for your study/project
% This example is from LaBGAS proj_erythritol_4a
% (https://gin.g-node.org/labgas/proj_erythritol_4a)
%
%
% DEPENDENCIES
%
% LaBGAScore Github repo on Matlab path, with subfolders
% https://github.com/labgas/LaBGAScore
%
%
% INPUTS
%
% Presentation .log files in sourcedata dir for each subject
%
%
% OUTPUTS
%
% events.tsv files for each run in BIDS dir for each subject
% phenotype.tsv file in BIDS/phenotype dir (optional)
%
%__________________________________________________________________________
%
% author: Lukas Van Oudenhove
% date:   December, 2021
%
%__________________________________________________________________________
% @(#)% LaBGAScore_prep_s1_write_events_tsv.m         v1.2        
% last modified: 2022/04/21


%% DEFINE DIRECTORIES, SUBJECTS, RUNS, CONDITIONS, AND IMPORT OPTIONS
%--------------------------------------------------------------------------

LaBGAScore_prep_s0_define_directories; % lukasvo edited from original LaBGAScore script to enable standalone functioning of proj_ery_4a dataset

subjs2write = {}; % enter subjects separated by comma if you only want to write files for selected subjects e.g. {'sub-01','sub-02'}
pheno_tsv = true; % turn to false if you do not wish to generate a phenotype.tsv file with trial-by-trial ratings; will only work if subjs2write is empty (i.e. when you loop over all your subjects)
pheno_name = 'ratings_online.tsv';

runnames = {'run-1','run-2','run-3','run-4','run-5','run-6'};
logfilenames = {'*_run1.log','*_run2.log','*_run3.log','*_run4.log','*_run5.log','*_run6.log'};
taskname = 'sweettaste_';
sweet_labels = {'sucrose delivery';'erythritol delivery';'sucralose delivery';'control delivery'}; % labels for sweet substance delivery in Code var of logfile
swallow_rinse_labels = {'sucrose_swallowing';'erythritol_swallowing';'sucralose_swallowing';'control_swallowing'}; % labels for swallowing cue presentation after sweet substance delivery in Code var of logfile
rating_labels = {'Sucrose','Erythritol','Sucralose','Control'}; % labels for start of rating period in Code var of logfile
fixation_labels = {'fixation_cross','Sucrose fixation cross','Erythritol fixation cross','Sucralose fixation cross','Control fixation cross'}; % labels for fixation cross in Code var of logfile
events_interest = {'sucrose','erythritol','sucralose','water'}; % names of events of interest to be written to events.tsv
events_nuisance = {'swallow_rinse','rating'}; % names of nuisance events to be written to events.tsv

varNames = {'Trial','Event Type','Code','Time','TTime','Uncertainty','Duration','Uncertainty','ReqTime','ReqDur','Stim Type','Pair Index'}; % varnames of logfile
selectedVarNames = [1:5 7 9]; % varnames we want to use in the script
varTypes = {'double','categorical','categorical','double','double','double','double','double','double','char','char','double'}; % matlab vartypes to be used when importing log file as table
delimiter = '\t';
dataStartLine = 5; % line on which actual data starts in logfile
extraColRule = 'ignore';

opts = delimitedTextImportOptions('VariableNames',varNames,...
                                'SelectedVariableNames',selectedVarNames,...
                                'VariableTypes',varTypes,...
                                'Delimiter',delimiter,...
                                'DataLines', dataStartLine,...
                                'ExtraColumnsRule',extraColRule); 


%% LOOP OVER SUBJECTS TO READ LOGFILES, CREATE TABLE WITH ONSETS AND DURATIONS, AND SAVE AS EVENTS.TSV TO BIDSSUBJDIRS
%-------------------------------------------------------------------------------------------------------------------
if ~isempty(subjs2write)
    [C,ia,~] = intersect(sourcesubjs,subjs2write);
    
    if ~isequal(C',subjs2write)
        error('\nsubject %s present in subjs2smooth not present in %s, please check before proceeding',subjs2smooth{~ismember(subjs2smooth,C)},derivdir);
    
    else
        
        for sub = ia'
        
        % DEFINE SUBJECT LEVEL DIRS
        subjsourcedir = sourcesubjdirs{sub};
        subjBIDSdir = fullfile(BIDSsubjdirs{sub},'func');

            % LOOP OVER RUNS
            for run = 1:size(logfilenames,2)
                
                logfilename = dir(fullfile(subjsourcedir,'logfiles',logfilenames{run}));
                logfilename = char(logfilename(:).name);
                logfilepath = fullfile(subjsourcedir,'logfiles',logfilename);
                
                if ~isfile(logfilepath)
                    warning('\nlogfile missing for run %d in %s, please check before proceeding',run,logfilepath);
                    continue
                
                elseif size(logfilepath,1) > 1
                    error('\nmore than one logfile with run index %s for %s, please check before proceeding',run,sourcesubjs{sub})
                
                else
                    log = readtable(logfilepath,opts);
                    log = log(~isnan(log.Trial),:);
                    time_zero = log.Time(log.Trial == 0 & log.EventType == 'Pulse'); % time for onsets and durations is counted from the first scanner pulse onwards
                        
                        if size(time_zero,1) > 1
                            error('\nambiguity about time zero in %s%s, please check logfile',subjs{sub},logfilenames{run});
                        end
                        
                    log.TimeZero = log.Time - time_zero;
                    log.onset = log.TimeZero ./ 10000; % convert to seconds
                    log(log.EventType == 'Pulse',:) = [];
                    log.trial_type = cell(height(log),1);
                        
                        for k = 1:height(log)
                            if ismember(log.Code(k),swallow_rinse_labels)
                                log.trial_type{k} = events_nuisance{1};
                            
                            elseif ismember(log.Code(k),rating_labels)
                                log.trial_type{k} = events_nuisance{2};
                                log.onset(k) = log.onset(k)+4;
                            
                            elseif ismember(log.Code(k),sweet_labels)
                                idx = (log.Code(k) == sweet_labels);
                                log.trial_type{k} = events_interest{idx'};
                            
                            elseif ismember(log.Code(k),fixation_labels)
                                log.trial_type{k} = 'fixation';
                            
                            else
                                log.trial_type{k} = '';
                            
                            end
                        end
                        
                    log.rating = zeros(height(log),1);
                    
                        for l = 1:height(log)
                            
                            if contains(char(log.Code(l)),'score','IgnoreCase',true)
                                scorestring = char(log.Code(l));
                                
                                if strcmp(scorestring(1,end-3:end),'-100')
                                    log.rating(l) = str2double(scorestring(1,end-3:end));
                                
                                elseif ~contains(scorestring(1,end-2:end),':')
                                    log.rating(l) = str2double(strtrim(scorestring(1,end-2:end)));
                                
                                else
                                    log.rating(l) = str2double(scorestring(1,end));
                                    
                                end
                            
                            else
                                log.rating(l) = NaN;
                                
                            end
                            
                        end
                        
                    log.trial_type = categorical(log.trial_type);
                    
                    log = log((~isundefined(log.trial_type) | ~isnan(log.rating)),:);
                    
                        for n = 1:height(log)
                            
                            if ismember(log.trial_type(n),events_interest)
                                log.rating(n) = log.rating(n+3);
                            
                            end
                            
                        end
                          
                    log = removevars(log,{'Trial','EventType','Code','Time','TTime','Duration','ReqTime','TimeZero'}); % get rid of junk variables from logfile we don't need
                    
                    log = log(~isundefined(log.trial_type),:);  
                    
                    log.duration = zeros(height(log),1);
                        
                        for m = 1:height(log)
                            
                            if ~isequal(log.trial_type(m),'fixation')
                                log.duration(m) = log.onset(m+1) - log.onset(m);
                            
                            else
                                log.duration(m) = NaN;
                            
                            end
                        end
                        
                    log = log(~isnan(log.duration),:);
                        
                    filename = fullfile(subjBIDSdir,[sourcesubjs{sub},'_task-',taskname,runnames{run},'_events.tsv']);
                    writetable(log,filename,'Filetype','text','Delimiter','\t');
                    clear logfile log time_zero filename

                end % if loop checking whether logfile exists

            end % for loop runs
            
        end % for loop subjects
    
    end % if loop checking subjs2write present in sourcesubjs        
    
else
    
    if ~isequal(sourcesubjs, BIDSsubjs)
        [D,~,~] = intersect(sourcesubjs,BIDSsubjs);
        error('\nsubject %s present in %s not present in %s, please check before proceeding',BIDSsubjs{~ismember(BIDSsubjs,D)},BIDSdir,derivdir);
        
    else
        
        if pheno_tsv
            pheno_file = table();
            pheno_dir = fullfile(BIDSdir,'phenotype');
                if ~isfolder(pheno_dir)
                    mkdir(pheno_dir);
                end
        end
        
        for sub = 1:size(sourcesubjs,1)
            
            if pheno_tsv
                pheno_file_subj = table();
            end

        % DEFINE SUBJECT LEVEL DIRS
        subjsourcedir = sourcesubjdirs{sub};
        subjBIDSdir = fullfile(BIDSsubjdirs{sub},'func');

            % LOOP OVER RUNS
            for run = 1:size(logfilenames,2)
                
                logfilename = dir(fullfile(subjsourcedir,'logfiles',logfilenames{run}));
                logfilename = char(logfilename(:).name);
                logfilepath = fullfile(subjsourcedir,'logfiles',logfilename);
                
                if ~isfile(logfilepath)
                    warning('\nlogfile missing for run %d in %s, please check before proceeding',run,logfilepath);
                    continue
                
                elseif size(logfilepath,1) > 1
                    error('\nmore than one logfile with run index %s for %s, please check before proceeding',run,sourcesubjs{sub})
                
                else
                    log = readtable(logfilepath,opts);
                    log = log(~isnan(log.Trial),:);
                    time_zero = log.Time(log.Trial == 0 & log.EventType == 'Pulse'); % time for onsets and durations is counted from the first scanner pulse onwards
                        
                        if size(time_zero,1) > 1
                            error('\nambiguity about time zero in %s%s, please check logfile',subjs{sub},logfilenames{run});
                        end
                        
                    log.TimeZero = log.Time - time_zero;
                    log.onset = log.TimeZero ./ 10000; % convert to seconds 
                    log(log.EventType == 'Pulse',:) = [];
                    log.trial_type = cell(height(log),1);
                        
                        for k = 1:height(log)
                            
                            if ismember(log.Code(k),swallow_rinse_labels)
                                log.trial_type{k} = events_nuisance{1};
                            
                            elseif ismember(log.Code(k),rating_labels)
                                log.trial_type{k} = events_nuisance{2};
                                log.onset(k) = log.onset(k)+4;
                            
                            elseif ismember(log.Code(k),sweet_labels)
                                idx = (log.Code(k) == sweet_labels);
                                log.trial_type{k} = events_interest{idx'};
                            
                            elseif ismember(log.Code(k),fixation_labels)
                                log.trial_type{k} = 'fixation';
                            
                            else
                                log.trial_type{k} = '';
                            
                            end
                        end
                    
                    log.rating = zeros(height(log),1);
                    
                        for l = 1:height(log)
                            
                            if contains(char(log.Code(l)),'score','IgnoreCase',true)
                                scorestring = char(log.Code(l));
                                
                                if strcmp(scorestring(1,end-3:end),'-100')
                                    log.rating(l) = str2double(scorestring(1,end-3:end));
                                
                                elseif ~contains(scorestring(1,end-2:end),':')
                                    log.rating(l) = str2double(strtrim(scorestring(1,end-2:end)));
                                
                                else
                                    log.rating(l) = str2double(scorestring(1,end));
                                    
                                end
                            
                            else
                                log.rating(l) = NaN;
                                
                            end
                            
                        end
                        
                    log.trial_type = categorical(log.trial_type);
                    
                    log = log((~isundefined(log.trial_type) | ~isnan(log.rating)),:);
                    
                        for n = 1:height(log)
                            
                            if ismember(log.trial_type(n),events_interest)
                                log.rating(n) = log.rating(n+3);
                            
                            end
                            
                        end
                          
                    log = removevars(log,{'Trial','EventType','Code','Time','TTime','Duration','ReqTime','TimeZero'}); % get rid of junk variables from logfile we don't need
                    
                    log = log(~isundefined(log.trial_type),:);  
                    
                    log.duration = zeros(height(log),1);
                        
                        for m = 1:height(log)
                            
                            if ~isequal(log.trial_type(m),'fixation')
                                log.duration(m) = log.onset(m+1) - log.onset(m);
                            
                            else
                                log.duration(m) = NaN;
                            
                            end
                        end
                        
                    log = log(~isnan(log.duration),:);
                        
                    filename = fullfile(subjBIDSdir,[sourcesubjs{sub},'_task-',taskname,runnames{run},'_events.tsv']);
                    writetable(log,filename,'Filetype','text','Delimiter','\t');
                    
                        if pheno_tsv
                            
                            log2 = log(~isnan(log.rating),:);
                            log2 = removevars(log2,{'onset','duration'});

                            for n = 1:height(log2)
                                log2.participant_id(n,:) = BIDSsubjs{sub};
                                log2.run_id(n) = run;
                                log2.trial_id_run(n) = n; % generates consecutive trial numbers within each run
                                log2.trial_id_concat(n) = height(pheno_file_subj) + n; % generates consecutive trial numbers over all conditions & runs
                            end  
                            
                            for o = 1:size(events_interest,2)
                                idx_run = log2.trial_type == events_interest{o};
                                log2.trial_id_cond_run(idx_run) = 1:sum(idx_run);
                                    if height(pheno_file_subj) > 0
                                        idx_sub = pheno_file_subj.trial_type == events_interest{o};
                                        log2.trial_id_cond_concat(idx_run) = sum(idx_sub)+1:(sum(idx_sub)+sum(idx_run));
                                    else
                                        log2.trial_id_cond_concat = log2.trial_id_cond_run;
                                    end
                                    clear idx_run idx_sub
                            end
                            
                            pheno_file_subj = [pheno_file_subj;log2];

                        end
                    
                    clear logfile log time_zero filename log2

                end % if loop checking whether logfile exists

            end % for loop runs
            
            pheno_file = [pheno_file;pheno_file_subj];
            clear pheno_file_subj;

        end % for loop subjects
        
        pheno_filename = fullfile(pheno_dir,pheno_name);
        writetable(pheno_file,pheno_filename,'Filetype','text','Delimiter','\t');
    
    end % if loop checking sourcesubjs == BIDSsubjs
    
end % if loop checking writing option