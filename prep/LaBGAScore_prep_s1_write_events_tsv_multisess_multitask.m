%% LaBGAScore_prep_s1_write_events_tsv_multisess_multitask
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
% This example is from LaBGAS proj_bitter-reward
% (https://gin.g-node.org/labgas/proj_bitter-reward)
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
% date:   January, 2023
%
%__________________________________________________________________________
% @(#)% LaBGAScore_prep_s1_write_events_tsv_multisess_multitask.m    v1.1        
% last modified: 2025/07/02
%
% adapted from LCN12_JULIE_first_level_analysis script by Patrick Dupont
%


%% DEFINE DIRECTORIES, SUBJECTS, RUNS, CONDITIONS, AND IMPORT OPTIONS
%--------------------------------------------------------------------------

study_prefix = ''; % STUDY-SPECIFIC

eval([study_prefix '_prep_s0_define_directories']);

subjs2write = {}; % enter subjects separated by comma if you only want to write files for selected subjects e.g. {'sub-01','sub-02'}
pheno_tsv = true; % turn to false if you do not wish to generate a phenotype.tsv file with trial-by-trial ratings; will only work if subjs2write is empty (i.e. when you loop over all your subjects)
pheno_name = 'food_images_wanting_liking.tsv';
nr_sess = 2;

runnames = {'run-1','run-2','run-1','run-2'}; % for naming of run in events.tsv files - two runs for each task
logfilenames = {'*_run1_photo_results.txt','*_run2_photo_results.txt','*_run3-FID_1_final*.log','*_run4-FID_1_final*.log'};
taskname1 = 'food_images_';
taskname2 = 'FID_';
image_labels = {'High','Low','Neutral','flow'}; % labels for sweet substance delivery in Code var of logfile
image_trial_types = {'rest','high calorie','low calorie','neutral','rating'};
image_trial_types_interest = {'high calorie','low calorie','neutral'};
FID_labels = {'open_circle','one_line','two_lines','feedback_10SP','feedback_2SP','feedback_0SP','right_arrow','left_arrow','11','12','99'};
FID_trial_types = {'cue C0','cue C2','cue C10','feedback_C10_win','feedback_C10_nowin','feedback_C2_win','feedback_C2_nowin','feedback_C0'};

vasfilename = '*_ratings_results.txt';

varNames = {'version','timestamp','block number','trial number','lost msec','free VRAM','rating','conditie','time_start_VAS','time_end_VAS','trial contents','unpleasantness','hunger','pfc','nausea','satiety','fullness','pl_eating','pl_pic','taste_pic'}; % varnames of logfile
selectedVarNames = [3:4 7:10 19:20]; % varnames we want to use in the script
varTypes = {'double','char','double','double','char','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double'}; % matlab vartypes to be used when importing log file as table
delimiter = '\t';
dataStartLine = 2; % line on which actual data starts in logfile
extraColRule = 'ignore';

opts = delimitedTextImportOptions('VariableNames',varNames,...
    'SelectedVariableNames',selectedVarNames,...
    'VariableTypes',varTypes,...
    'Delimiter',delimiter,...
    'DataLines', dataStartLine,...
    'ExtraColumnsRule',extraColRule); 

% 1 = feedback_C10_win, 2 = feedback_C10_nowin, 3 = feedback_C2_win
% 4 = feedback_C2_nowin, 5 = feedback_C0
order_trials = [5 4 3 3 1 4 1 5 2 5 ...
                3 3 1 3 2 5 3 2 5 1 ... 
                3 3 1 4 4 5 3 1 3 5 ...
                4 1 5 4 3 5 5 1 1 2 ...
                1 5 2 5 3 2 4 5 5 1 ...
                2 2 5 1 3];

            
%% LOOP OVER SUBJECTS TO READ LOGFILES, CREATE TABLE WITH ONSETS AND DURATIONS, AND SAVE AS EVENTS.TSV TO BIDSSUBJDIRS
%-------------------------------------------------------------------------------------------------------------------

if ~isempty(subjs2write)
    [C,ia,~] = intersect(sourcesubjs,subjs2write);
    
    if ~isequal(C',subjs2write)
        error('\nsubject %s present in subjs2write not present in %s, please check before proceeding',subjs2write{~ismember(subjs2write,C)},derivdir);
    
    else
        
        % LOOP OVER SUBJECTS
        
        for sub = ia'
        
            % LOOP OVER SESSIONS
            
            for ses = 1:nr_sess
        
            % DEFINE SUBJECT & SESSION LEVEL DIRS
            
            subjsourcedir = sourcesubjdirs{sub};
            [~,subjectcode,~] = fileparts(subjsourcedir);
            subjBIDSdir = fullfile(BIDSsubjdirs{sub});
            sessubjsourcedir = fullfile(subjsourcedir,['ses-' num2str(ses)]);
            sessubjBIDSdir = fullfile(subjBIDSdir,['ses-' num2str(ses)],'func');
            logsessubjsourcedir = fullfile(sessubjsourcedir,'logfiles');
            
                % LOOP OVER RUNS
                
                for run = 1:size(runnames,2)
                    
                    % DEFINE LOGFILES

                    logfilename = dir(fullfile(logsessubjsourcedir,logfilenames{run}));
                    logfilename = char(logfilename(:).name);
                    logfilepath = fullfile(logsessubjsourcedir,logfilename);
                    
                    % SANITY CHECK

                    if ~isfile(logfilepath)
                        warning('\nlogfile missing for run %d in %s, please check before proceeding',run,logfilepath);
                        continue

                    elseif size(logfilepath,1) > 1
                        error('\nmore than one logfile with run index %s for %s, please check before proceeding',run,sourcesubjs{sub})

                    else
                        
                        % FOOD IMAGES RUNS
                        
                        if run < 3 
                            
                            % read onsets file - Patrick's code
                            
                            fid = fopen(logfilepath,'r');
                            datafile = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %*[^\n]','delimiter','\t');
                            fclose(fid);
                            
                            % first row of datafile is a header and should be neglected
                            % column 8 of datafile is the start time in ms since start of fMRI scan
                            % column 9 of datafile is the end time in ms since start of fMRI scan
                            % column 10 of datafile is the image

                            tmp8 = datafile{8};
                            tmp8 = tmp8(2:end);
                            % make tmp8 an array of onsets
                            onsets = zeros(length(tmp8),1);
                            
                                for i = 1:length(tmp8)
                                   onsets(i) = str2num(tmp8{i});
                                end
                               
                            tmp10 = datafile{10};
                            tmp10 = tmp10(2:end); % this is the list of images
                            index_LC     = [];
                            index_HC     = [];
                            index_NF     = [];
                            index_rest   = [];
                            index_rating = [];
                            
                               for i = 1:length(tmp10)
                                   
                                  if strfind(tmp10{i},image_labels{1}) == 1
                                      index_HC = [index_HC i];
                                   elseif strfind(tmp10{i},image_labels{2}) == 1
                                      index_LC = [index_LC i];
                                   elseif strfind(tmp10{i},image_labels{3}) == 1
                                      index_NF = [index_NF i];
                                   elseif strfind(tmp10{i},image_labels{4}) == 1
                                      index_rest = [index_rest i];
                                  end
                                   
                               end

                            % since 8 images were shown, the index_LC, index_HC and index_NF
                            % contain always blocks of 8 subsequent images. We only need the
                            % onset of the block of images.

                            onsets_blocks{1} = onsets(index_rest)/1000; % in seconds
                            onsets_blocks{2} = onsets(index_HC(1:8:end))/1000; % in seconds
                            onsets_blocks{3} = onsets(index_LC(1:8:end))/1000; % in seconds
                            onsets_blocks{4} = onsets(index_NF(1:8:end))/1000; % in seconds
                           
                            tmp1 = onsets(index_HC(8:8:end))/1000 + 3; % onset of the rating block after HC
                            tmp1_durations = 12*ones(length(tmp1),1);
                            tmp2 = onsets(index_LC(8:8:end))/1000 + 3; % onset of the rating block after LC 
                            tmp2_durations = 12*ones(length(tmp2),1);
                            tmp3 = onsets(index_NF(8:8:end))/1000 + 3; % onset of the rating block after NF
                            tmp3_durations = 6*ones(length(tmp3),1);
                            tmp = [tmp1; tmp2; tmp3];
                            tmp_durations = [tmp1_durations; tmp2_durations; tmp3_durations];
                            [onsets_blocks{5},listorder] = sort(tmp,'ascend');
                           
                            durations{1} = 24; % rest period of 24s
                            durations{2} = 24; % 8 HC images x 3s
                            durations{3} = 24; % 8 LC images x 3s
                            durations{4} = 24; % 8 NF images x 3s        
                            durations{5} = tmp_durations(listorder);
                           
                            % conversion to events.tsv table format - Lukas' code
                            
                            onset = [];
                            duration = [];
                            trial_type = categorical();
                            
                                for cond = 1:size(image_trial_types,2)
                                    
                                    if cond < 5
                                        trial_type_block = categorical();
                                        duration_block = [];
                                        
                                        for block = 1:size(onsets_blocks{cond},1)
                                            trial_type_block(block,1) = image_trial_types{cond};
                                            duration_block(block,1) = durations{cond};
                                        end
                                        
                                        onset = [onset; onsets_blocks{cond}];
                                        duration = [duration; duration_block];
                                        trial_type = [trial_type; trial_type_block];
                                        clear trial_type_block duration_block
                                        
                                    else
                                        trial_type_block = categorical();
                                        
                                        for block = 1:size(onsets_blocks{cond},1)
                                            trial_type_block(block,1) = image_trial_types{cond};
                                        end
                                        
                                        onset = [onset; onsets_blocks{cond}];
                                        duration = [duration; durations{cond}];
                                        trial_type = [trial_type; trial_type_block];
                                        clear trial_type_block
                                        
                                    end % if loop checking condition type
                                    
                                end % for loop over conditions
                                
                            log=table(trial_type,onset,duration);
                            log=sortrows(log,'onset','ascend');
                            
                            % read ratings file as a table - Lukas' code
                            
                            ratingfilename = dir(fullfile(logsessubjsourcedir,vasfilename));
                            ratingfilename = char(ratingfilename(run).name);
                            ratingfilepath = fullfile(logsessubjsourcedir,ratingfilename);

                                if ~isfile(ratingfilepath)
                                    warning('\nlogfile missing for run %d in %s, please check before proceeding',run,ratingfilepath);
                                    continue

                                elseif size(ratingfilepath,1) > 1
                                    error('\nmore than one logfile with run index %s for %s, please check before proceeding',run,sourcesubjs{sub})

                                else
                            
                                    vas=readtable(ratingfilepath,opts);
                                    vas.wanting = (vas.rating).*(vas.pl_pic);
                                    vas.liking = (vas.rating).*(vas.taste_pic);
                                    vas.keep = vas.pl_pic + vas.taste_pic;
                                    vas = vas((vas.keep == 1),:);
                                    idx_truezero = (vas.rating == 0);
                                    
                                    wanting = vas.wanting(1:end-1);
                                    liking = vas.liking(2:end); % wanting always before liking, we want them on the same row, as well as the index
                                    ratings = table(wanting,liking);    
                                        if sum(idx_truezero) > 0
                                            if vas.pl_pic(idx_truezero) == 1
                                                idx_truezero = idx_truezero(1:end-1);
                                            elseif vas.taste_pic(idx_truezero) == 1
                                                idx_truezero = idx_truezero(2:end);
                                            end
                                            idx_zero = logical((ratings.wanting + ratings.liking > 0) + idx_truezero);
                                        else
                                            idx_zero = logical(ratings.wanting + ratings.liking > 0);
                                        end
                                    ratings = ratings(idx_zero,:);
                                    idx_rating = ismember(log.trial_type,image_trial_types_interest);
                                    
                                    log.wanting(idx_rating) = ratings.wanting;
                                    log.liking(idx_rating) = ratings.liking;
                                    log.wanting(log.wanting == 0) = NaN;
                                    log.liking(log.liking == 0) = NaN;
                                    
                                end
                            
                            % write events.tsv file - Lukas' code
                            
                            filename = fullfile(sessubjBIDSdir,[sourcesubjs{sub},['_ses-' num2str(ses)],'_task-',taskname1,runnames{run},'_events.tsv']);
                            writetable(log,filename,'Filetype','text','Delimiter','\t');
                            clear datafile index_* tmp* log onset* duration* vas wanting liking ratings idx_rating trial_type filename
                            
                        else % FID RUNS
                            
                            % read onsets file - Patrick's code
                            
                            % read headers
                            
                            fid = fopen(logfilepath,'r'); 
                            tmp = textscan(fid,'%s',1,'delimiter','\n'); 
                            creationtime = textscan(fid,'%s',1,'delimiter','\n'); 
                            tmp = textscan(fid,'%s',1,'delimiter','\n'); 
                            headerstring = textscan(fid,'%s',13,'delimiter','\t'); 
                           
                            % headerstring should now contain the header of the log file
                            
                            %    'Subject'
                            %    'Trial'
                            %    'Event Type'
                            %    'Code'
                            %    'Time'
                            %    'TTime'
                            %    'Uncertainty'
                            %    'Duration'
                            %    'Uncertainty'
                            %    'ReqTime'
                            %    'ReqDur'
                            %    'Stim Type'
                            %    'Pair Index'
                           
                            tmp = textscan(fid,'%s',1,'delimiter','\n');  % empty line        
                           
                            % read now the data part
                            
                            datafile = textscan(fid,'%s%d%s%s%f%*[^\n]','delimiter','\t');
                            fclose(fid);
                            
                            %    condition          code in logfile
                            %    'cue C0'           open_circle
                            %    'cue C2'           one_line
                            %    'cue C10'          two_lines
                            %    'cor after C10'    feedback_10SP AND [(right_arrow followed by response 12) OR (left_arrow followed by response 11)]
                            %    'incor after C10'  feedback_10SP AND [(right_arrow followed by response 12) OR (left_arrow followed by response 11)]
                            %    'cor after C2'     feedback_2SP AND [(right_arrow followed by response 12) OR (left_arrow followed by response 11)] 
                            %    'incor after C2'   feedback_2SP AND [(right_arrow followed by response 11) OR (left_arrow followed by response 12)]
                            %    'after C0'         feedback_0SP AND [(right_arrow followed by response 11 or 12) OR (left_arrow followed by response 11 or 12)]

                            % get the onsets for each condition
                            
                            timing = datafile{5}/10000; % in seconds
                            code = datafile{4};
                            index_cue_C0     = [];
                            index_cue_C2     = [];
                            index_cue_C10    = [];
                            index_feedback_C10 = [];
                            index_feedback_C2  = [];
                            index_feedback_C0  = [];
                            index_Rarrow  = [];
                            index_Larrow  = [];
                            index_responseL = [];
                            index_responseR = [];
                            index_scan = []; % all the pulses from the scanner
                           
                            for i = 1:length(code)
                               if strcmp(code{i},FID_labels{1})
                                  index_cue_C0 = [index_cue_C0 i]; 
                               elseif strcmp(code{i},FID_labels{2})
                                  index_cue_C2 = [index_cue_C2 i]; 
                               elseif strcmp(code{i},FID_labels{3})
                                  index_cue_C10 = [index_cue_C10 i]; 
                               elseif strcmp(code{i},FID_labels{4})
                                  index_feedback_C10 = [index_feedback_C10 i];
                               elseif strcmp(code{i},FID_labels{5})
                                  index_feedback_C2 = [index_feedback_C2 i];
                               elseif strcmp(code{i},FID_labels{6})
                                  index_feedback_C0 = [index_feedback_C0 i];
                               elseif strcmp(code{i},FID_labels{7})
                                  index_Rarrow = [index_Rarrow i];
                               elseif strcmp(code{i},FID_labels{8})
                                  index_Larrow = [index_Larrow i];
                               elseif strcmp(code{i},FID_labels{9})
                                  index_responseL = [index_responseL i];
                               elseif strcmp(code{i},FID_labels{10})
                                  index_responseR = [index_responseR i];
                               elseif strcmp(code{i},FID_labels{11})
                                  index_scan = [index_scan i];
                               end
                            end
                           
                            index_trials_start_cue = sort([index_cue_C0 index_cue_C2 index_cue_C10]); 
                           
                           % determine the onsets taken into account that the logfile
                           % starts when the dummy scan started. So we have to correct the
                           % timing with respect to the onset of the fifth scan.
                           
                           clear onsets
                           onsets{1} = timing(index_cue_C0)-timing(index_scan(5)); 
                           onsets{2} = timing(index_cue_C2)-timing(index_scan(5));
                           onsets{3} = timing(index_cue_C10)-timing(index_scan(5));
                           
                           timing_feedback_0SP = timing(index_feedback_C0)-timing(index_scan(5));
                           timing_feedback_2SP = timing(index_feedback_C2)-timing(index_scan(5));
                           timing_feedback_10SP = timing(index_feedback_C10)-timing(index_scan(5));
                           timing_Larrow = timing(index_Larrow)-timing(index_scan(5));
                           timing_Rarrow = timing(index_Rarrow)-timing(index_scan(5));
                           timing_Lresponse = timing(index_responseL)-timing(index_scan(5));
                           timing_Rresponse = timing(index_responseR)-timing(index_scan(5));
                           timing_feedback = timing(index_trials_start_cue)-timing(index_scan(5))+0.750+3+1; % 0.750s cue display, 3s delay, 1s target

                           % determine which one corresponds to a correct or incorrect response
                           
                           onsets{4} = []; % 'cor after C10' feedback_10SP AND [(right_arrow followed by response 12) OR (left_arrow followed by response 11)]
                           onsets{5} = []; % 'incor after C10' feedback_10SP AND [(right_arrow followed by response 12) OR (left_arrow followed by response 11)]
                           onsets{6} = []; % 'cor after C2' feedback_2SP AND [(right_arrow followed by response 12) OR (left_arrow followed by response 11)]
                           onsets{7} = []; % 'incor after C2' feedback_2SP AND [(right_arrow followed by response 11) OR (left_arrow followed by response 12)]
                           onsets{8} = []; % 'after C0' feedback_0SP AND [(right_arrow followed by response 11 or 12) OR (left_arrow followed by response 11 or 12)]
                           
                               for i = 1:length(timing_feedback_10SP)
                                   time0 = timing_feedback_10SP(i);
                                   La = intersect(find(timing_Larrow < time0),find(timing_Larrow > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
                                   Ra = intersect(find(timing_Rarrow < time0),find(timing_Rarrow > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
                                   Lr = intersect(find(timing_Lresponse < time0),find(timing_Lresponse > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
                                   Rr = intersect(find(timing_Rresponse < time0),find(timing_Rresponse > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
                                   
                                   if (length(La) == 1 && length(Lr) == 1) || (length(Ra) == 1 && length(Rr) == 1) % correct response
                                      onsets{4} = [onsets{4}; time0];
                                   end
                                   
                               end
                               
                               for i = 1:length(timing_feedback_2SP)
                                   time0 = timing_feedback_2SP(i);
                                   La = intersect(find(timing_Larrow < time0),find(timing_Larrow > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
                                   Ra = intersect(find(timing_Rarrow < time0),find(timing_Rarrow > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
                                   Lr = intersect(find(timing_Lresponse < time0),find(timing_Lresponse > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
                                   Rr = intersect(find(timing_Rresponse < time0),find(timing_Rresponse > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
                                   
                                   if (length(La) == 1 && length(Lr) == 1) || (length(Ra) == 1 && length(Rr) == 1) % correct response
                                      % now check if this is a win trial or a no win trial
                                      onsets{6} = [onsets{6}; time0];
                                   end
                                   
                               end
                               
                               for j = 1:length(timing_feedback_0SP)
                                   time0 = timing_feedback_0SP(j);
                                   La = intersect(find(timing_Larrow < time0),find(timing_Larrow > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
                                   Ra = intersect(find(timing_Rarrow < time0),find(timing_Rarrow > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
                                   Lr = intersect(find(timing_Lresponse < time0),find(timing_Lresponse > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
                                   Rr = intersect(find(timing_Rresponse < time0),find(timing_Rresponse > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
                                   if (length(La) == 1 && length(Lr) == 1) || (length(Ra) == 1 && length(Rr) == 1) % correct response
                                      % now check if this is a win trial or a no win trial
                                      index_trial = find(abs(timing_feedback-time0) == min(abs(timing_feedback-time0)));
                                      if order_trials(index_trial) == 2
                                         onsets{5} = [onsets{5}; time0];
                                      elseif order_trials(index_trial) == 4
                                         onsets{7} = [onsets{7}; time0];
                                      elseif order_trials(index_trial) == 5
                                         onsets{8} = [onsets{8}; time0];
                                      else
                                         fprintf('Something is wrong: inconsistencies found for feedback_C0 \n');
                                      end
                                   elseif (length(La) == 1 && length(Rr) == 1) || (length(Ra) == 1 && length(Lr) == 1)  % incorrect response
                                      % now check if this is a win trial or a no win trial
                                      index_trial = find(abs(timing_feedback-time0) == min(abs(timing_feedback-time0)));
                                      if order_trials(index_trial) == 2
                                         onsets{5} = [onsets{5}; time0];
                                      elseif order_trials(index_trial) == 4
                                         onsets{7} = [onsets{7}; time0];
                                      elseif order_trials(index_trial) == 5
                                         onsets{8} = [onsets{8}; time0];
                                      else
                                         fprintf('Something is wrong: inconsistencies found for feedback_C0 \n');
                                      end
                                      
                                   end
                                   
                               end
                           
                           % conversion to events.tsv table format - Lukas' code
                           
                           onset = [];
                           trial_type = categorical();
                            
                            for cond = 1:size(FID_trial_types,2)
                                
                                    trial_type_event = categorical();

                                    for ev = 1:size(onsets{cond},1)
                                        trial_type_event(ev,1) = FID_trial_types{cond};
                                    end

                                    onset = [onset; onsets{cond}];
                                    trial_type = [trial_type; trial_type_event];
                                    clear trial_type_event

                            end
                                
                           log=table(trial_type,onset);
                           log.duration = zeros(height(log),1);
                           log=sortrows(log,'onset','ascend');
                             

                           filename = fullfile(sessubjBIDSdir,[sourcesubjs{sub},['_ses-' num2str(ses)],'_task-',taskname2,runnames{run},'_events.tsv']);
                           writetable(log,filename,'Filetype','text','Delimiter','\t');
                           clear datafile index_* tmp* log onset* duration* trial_type filename
                           
                        end % if loop determining run number

                    end % if loop checking whether logfile exists

                end % for loop runs
            
            end % for loop sessions
            
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

            % LOOP OVER SESSIONS
            
            for ses = 1:2
                
            if pheno_tsv
                pheno_file_ses = table();
            end
        
            % DEFINE SUBJECT & SESSION LEVEL DIRS
            
            subjsourcedir = sourcesubjdirs{sub};
            [~,subjectcode,~] = fileparts(subjsourcedir);
            subjBIDSdir = fullfile(BIDSsubjdirs{sub});
            sessubjsourcedir = fullfile(subjsourcedir,['ses-' num2str(ses)]);
            sessubjBIDSdir = fullfile(subjBIDSdir,['ses-' num2str(ses)],'func');
            logsessubjsourcedir = fullfile(sessubjsourcedir,'logfiles');
            
                % LOOP OVER RUNS
                
                for run = 1:size(runnames,2)
                    
                    % DEFINE LOGFILES

                    logfilename = dir(fullfile(logsessubjsourcedir,logfilenames{run}));
                    logfilename = char(logfilename(:).name);
                    logfilepath = fullfile(logsessubjsourcedir,logfilename);
                    
                    % SANITY CHECK

                    if ~isfile(logfilepath)
                        warning('\nlogfile missing for run %d in %s, please check before proceeding',run,logfilepath);
                        continue

                    elseif size(logfilepath,1) > 1
                        error('\nmore than one logfile with run index %s for %s, please check before proceeding',run,sourcesubjs{sub})

                    else
                        
                        % FOOD IMAGES RUNS
                        
                        if run < 3 
                            
                            % read onsets file - Patrick's code
                            
                            fid = fopen(logfilepath,'r');
                            datafile = textscan(fid,'%s %s %s %s %s %s %s %s %s %s %*[^\n]','delimiter','\t');
                            fclose(fid);
                            
                            % first row of datafile is a header and should be neglected
                            % column 8 of datafile is the start time in ms since start of fMRI scan
                            % column 9 of datafile is the end time in ms since start of fMRI scan
                            % column 10 of datafile is the image

                            tmp8 = datafile{8};
                            tmp8 = tmp8(2:end);
                            % make tmp8 an array of onsets
                            onsets = zeros(length(tmp8),1);
                            
                                for i = 1:length(tmp8)
                                   onsets(i) = str2num(tmp8{i});
                                end
                               
                            tmp10 = datafile{10};
                            tmp10 = tmp10(2:end); % this is the list of images
                            index_LC     = [];
                            index_HC     = [];
                            index_NF     = [];
                            index_rest   = [];
                            index_rating = [];
                            
                               for i = 1:length(tmp10)
                                   
                                  if strfind(tmp10{i},image_labels{1}) == 1
                                      index_HC = [index_HC i];
                                   elseif strfind(tmp10{i},image_labels{2}) == 1
                                      index_LC = [index_LC i];
                                   elseif strfind(tmp10{i},image_labels{3}) == 1
                                      index_NF = [index_NF i];
                                   elseif strfind(tmp10{i},image_labels{4}) == 1
                                      index_rest = [index_rest i];
                                  end
                                   
                               end

                            % since 8 images were shown, the index_LC, index_HC and index_NF
                            % contain always blocks of 8 subsequent images. We only need the
                            % onset of the block of images.

                            onsets_blocks{1} = onsets(index_rest)/1000; % in seconds
                            onsets_blocks{2} = onsets(index_HC(1:8:end))/1000; % in seconds
                            onsets_blocks{3} = onsets(index_LC(1:8:end))/1000; % in seconds
                            onsets_blocks{4} = onsets(index_NF(1:8:end))/1000; % in seconds
                           
                            tmp1 = onsets(index_HC(8:8:end))/1000 + 3; % onset of the rating block after HC
                            tmp1_durations = 12*ones(length(tmp1),1);
                            tmp2 = onsets(index_LC(8:8:end))/1000 + 3; % onset of the rating block after LC 
                            tmp2_durations = 12*ones(length(tmp2),1);
                            tmp3 = onsets(index_NF(8:8:end))/1000 + 3; % onset of the rating block after NF
                            tmp3_durations = 6*ones(length(tmp3),1);
                            tmp = [tmp1; tmp2; tmp3];
                            tmp_durations = [tmp1_durations; tmp2_durations; tmp3_durations];
                            [onsets_blocks{5},listorder] = sort(tmp,'ascend');
                           
                            durations{1} = 24; % rest period of 24s
                            durations{2} = 24; % 8 HC images x 3s
                            durations{3} = 24; % 8 LC images x 3s
                            durations{4} = 24; % 8 NF images x 3s        
                            durations{5} = tmp_durations(listorder);
                           
                            % conversion to events.tsv table format - Lukas' code
                            
                            onset = [];
                            duration = [];
                            trial_type = categorical();
                            
                                for cond = 1:size(image_trial_types,2)
                                    
                                    if cond < 5
                                        trial_type_block = categorical();
                                        duration_block = [];
                                        
                                        for block = 1:size(onsets_blocks{cond},1)
                                            trial_type_block(block,1) = image_trial_types{cond};
                                            duration_block(block,1) = durations{cond};
                                        end
                                        
                                        onset = [onset; onsets_blocks{cond}];
                                        duration = [duration; duration_block];
                                        trial_type = [trial_type; trial_type_block];
                                        clear trial_type_block duration_block
                                        
                                    else
                                        trial_type_block = categorical();
                                        
                                        for block = 1:size(onsets_blocks{cond},1)
                                            trial_type_block(block,1) = image_trial_types{cond};
                                        end
                                        
                                        onset = [onset; onsets_blocks{cond}];
                                        duration = [duration; durations{cond}];
                                        trial_type = [trial_type; trial_type_block];
                                        clear trial_type_block
                                        
                                    end % if loop checking condition type
                                    
                                end % for loop over conditions
                                
                            log=table(trial_type,onset,duration);
                            log=sortrows(log,'onset','ascend');
                            
                            % read ratings file as a table - Lukas' code
                            
                            ratingfilename = dir(fullfile(logsessubjsourcedir,vasfilename));
                            ratingfilename = char(ratingfilename(run).name);
                            ratingfilepath = fullfile(logsessubjsourcedir,ratingfilename);

                                if ~isfile(ratingfilepath)
                                    warning('\nlogfile missing for run %d in %s, please check before proceeding',run,ratingfilepath);
                                    continue

                                elseif size(ratingfilepath,1) > 1
                                    error('\nmore than one logfile with run index %s for %s, please check before proceeding',run,sourcesubjs{sub})

                                else
                            
                                    vas=readtable(ratingfilepath,opts);
                                    vas.wanting = (vas.rating).*(vas.pl_pic);
                                    vas.liking = (vas.rating).*(vas.taste_pic);
                                    vas.keep = vas.pl_pic + vas.taste_pic;
                                    vas = vas((vas.keep == 1),:);
                                    idx_truezero = (vas.rating == 0);
                                    
                                    wanting = vas.wanting(1:end-1);
                                    liking = vas.liking(2:end); % wanting always before liking, we want them on the same row, as well as the index
                                    ratings = table(wanting,liking);    
                                        if sum(idx_truezero) > 0
                                            if vas.pl_pic(idx_truezero) == 1
                                                idx_truezero = idx_truezero(1:end-1);
                                            elseif vas.taste_pic(idx_truezero) == 1
                                                idx_truezero = idx_truezero(2:end);
                                            end
                                            idx_zero = logical((ratings.wanting + ratings.liking > 0) + idx_truezero);
                                        else
                                            idx_zero = logical(ratings.wanting + ratings.liking > 0);
                                        end
                                    ratings = ratings(idx_zero,:);
                                    idx_rating = ismember(log.trial_type,image_trial_types_interest);
                                    
                                    log.wanting(idx_rating) = ratings.wanting;
                                    log.liking(idx_rating) = ratings.liking;
                                    log.wanting(log.wanting == 0) = NaN;
                                    log.liking(log.liking == 0) = NaN;
                                    
                                end
                            
                            % write events.tsv file - Lukas' code
                            
                            filename = fullfile(sessubjBIDSdir,[sourcesubjs{sub},['_ses-' num2str(ses)],'_task-',taskname1,runnames{run},'_events.tsv']);
                            writetable(log,filename,'Filetype','text','Delimiter','\t');
                            clear datafile index_* tmp* onset* duration* vas wanting liking ratings idx_rating trial_type filename
                            
                            % write phenotype.tsv file if requested - Lukas' code
                            
                                if pheno_tsv

                                    log2 = log(~isnan(log.wanting),:);
                                    log2 = removevars(log2,{'onset','duration'});

                                    for n = 1:height(log2)
                                        log2.participant_id(n,:) = BIDSsubjs{sub};
                                        log2.session(n) = ses;
                                        log2.run_id(n) = run;
                                        log2.trial_id_run(n) = n; % generates consecutive trial numbers within each run
                                        log2.trial_id_concat(n) = height(pheno_file_subj) + n; % generates consecutive trial numbers over all conditions & runs
                                            if ses == 1
                                                log2.condition(n,:) = 'bitter';
                                            elseif ses == 2
                                                log2.condition(n,:) = 'saline';
                                            end
                                    end  

                                    for o = 1:size(image_trial_types_interest,2)
                                        idx_run = log2.trial_type == image_trial_types_interest{o};
                                        log2.trial_id_cond_run(idx_run) = 1:sum(idx_run);
                                            if height(pheno_file_subj) > 0
                                                idx_sub = pheno_file_subj.trial_type == image_trial_types_interest{o};
                                                log2.trial_id_cond_concat(idx_run) = sum(idx_sub)+1:(sum(idx_sub)+sum(idx_run));
                                            else
                                                log2.trial_id_cond_concat = log2.trial_id_cond_run;
                                            end
                                            clear idx_run idx_sub
                                    end

                                    pheno_file_ses = [pheno_file_ses;log2];

                                end % if loop checking whether pheno_tsv is requested
                                
                                clear datafile index_* tmp* log log2 onset* duration* vas wanting liking ratings idx_rating trial_type filename
                            
                        else % FID RUNS
                            
                            % read onsets file - Patrick's code
                            
                            % read headers
                            
                            fid = fopen(logfilepath,'r'); 
                            tmp = textscan(fid,'%s',1,'delimiter','\n'); 
                            creationtime = textscan(fid,'%s',1,'delimiter','\n'); 
                            tmp = textscan(fid,'%s',1,'delimiter','\n'); 
                            headerstring = textscan(fid,'%s',13,'delimiter','\t'); 
                           
                            % headerstring should now contain the header of the log file
                            
                            %    'Subject'
                            %    'Trial'
                            %    'Event Type'
                            %    'Code'
                            %    'Time'
                            %    'TTime'
                            %    'Uncertainty'
                            %    'Duration'
                            %    'Uncertainty'
                            %    'ReqTime'
                            %    'ReqDur'
                            %    'Stim Type'
                            %    'Pair Index'
                           
                            tmp = textscan(fid,'%s',1,'delimiter','\n');  % empty line        
                           
                            % read now the data part
                            
                            datafile = textscan(fid,'%s%d%s%s%f%*[^\n]','delimiter','\t');
                            fclose(fid);
                            
                            %    condition          code in logfile
                            %    'cue C0'           open_circle
                            %    'cue C2'           one_line
                            %    'cue C10'          two_lines
                            %    'cor after C10'    feedback_10SP AND [(right_arrow followed by response 12) OR (left_arrow followed by response 11)]
                            %    'incor after C10'  feedback_10SP AND [(right_arrow followed by response 12) OR (left_arrow followed by response 11)]
                            %    'cor after C2'     feedback_2SP AND [(right_arrow followed by response 12) OR (left_arrow followed by response 11)] 
                            %    'incor after C2'   feedback_2SP AND [(right_arrow followed by response 11) OR (left_arrow followed by response 12)]
                            %    'after C0'         feedback_0SP AND [(right_arrow followed by response 11 or 12) OR (left_arrow followed by response 11 or 12)]

                            % get the onsets for each condition
                            
                            timing = datafile{5}/10000; % in seconds
                            code = datafile{4};
                            index_cue_C0     = [];
                            index_cue_C2     = [];
                            index_cue_C10    = [];
                            index_feedback_C10 = [];
                            index_feedback_C2  = [];
                            index_feedback_C0  = [];
                            index_Rarrow  = [];
                            index_Larrow  = [];
                            index_responseL = [];
                            index_responseR = [];
                            index_scan = []; % all the pulses from the scanner
                           
                            for i = 1:length(code)
                               if strcmp(code{i},FID_labels{1})
                                  index_cue_C0 = [index_cue_C0 i]; 
                               elseif strcmp(code{i},FID_labels{2})
                                  index_cue_C2 = [index_cue_C2 i]; 
                               elseif strcmp(code{i},FID_labels{3})
                                  index_cue_C10 = [index_cue_C10 i]; 
                               elseif strcmp(code{i},FID_labels{4})
                                  index_feedback_C10 = [index_feedback_C10 i];
                               elseif strcmp(code{i},FID_labels{5})
                                  index_feedback_C2 = [index_feedback_C2 i];
                               elseif strcmp(code{i},FID_labels{6})
                                  index_feedback_C0 = [index_feedback_C0 i];
                               elseif strcmp(code{i},FID_labels{7})
                                  index_Rarrow = [index_Rarrow i];
                               elseif strcmp(code{i},FID_labels{8})
                                  index_Larrow = [index_Larrow i];
                               elseif strcmp(code{i},FID_labels{9})
                                  index_responseL = [index_responseL i];
                               elseif strcmp(code{i},FID_labels{10})
                                  index_responseR = [index_responseR i];
                               elseif strcmp(code{i},FID_labels{11})
                                  index_scan = [index_scan i];
                               end
                            end
                           
                            index_trials_start_cue = sort([index_cue_C0 index_cue_C2 index_cue_C10]); 
                           
                           % determine the onsets taken into account that the logfile
                           % starts when the dummy scan started. So we have to correct the
                           % timing with respect to the onset of the fifth scan.
                           
                           clear onsets
                           onsets{1} = timing(index_cue_C0)-timing(index_scan(5)); 
                           onsets{2} = timing(index_cue_C2)-timing(index_scan(5));
                           onsets{3} = timing(index_cue_C10)-timing(index_scan(5));
                           
                           timing_feedback_0SP = timing(index_feedback_C0)-timing(index_scan(5));
                           timing_feedback_2SP = timing(index_feedback_C2)-timing(index_scan(5));
                           timing_feedback_10SP = timing(index_feedback_C10)-timing(index_scan(5));
                           timing_Larrow = timing(index_Larrow)-timing(index_scan(5));
                           timing_Rarrow = timing(index_Rarrow)-timing(index_scan(5));
                           timing_Lresponse = timing(index_responseL)-timing(index_scan(5));
                           timing_Rresponse = timing(index_responseR)-timing(index_scan(5));
                           timing_feedback = timing(index_trials_start_cue)-timing(index_scan(5))+0.750+3+1; % 0.750s cue display, 3s delay, 1s target

                           % determine which one corresponds to a correct or incorrect response
                           
                           onsets{4} = []; % 'cor after C10' feedback_10SP AND [(right_arrow followed by response 12) OR (left_arrow followed by response 11)]
                           onsets{5} = []; % 'incor after C10' feedback_10SP AND [(right_arrow followed by response 12) OR (left_arrow followed by response 11)]
                           onsets{6} = []; % 'cor after C2' feedback_2SP AND [(right_arrow followed by response 12) OR (left_arrow followed by response 11)]
                           onsets{7} = []; % 'incor after C2' feedback_2SP AND [(right_arrow followed by response 11) OR (left_arrow followed by response 12)]
                           onsets{8} = []; % 'after C0' feedback_0SP AND [(right_arrow followed by response 11 or 12) OR (left_arrow followed by response 11 or 12)]
                           
                               for i = 1:length(timing_feedback_10SP)
                                   time0 = timing_feedback_10SP(i);
                                   La = intersect(find(timing_Larrow < time0),find(timing_Larrow > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
                                   Ra = intersect(find(timing_Rarrow < time0),find(timing_Rarrow > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
                                   Lr = intersect(find(timing_Lresponse < time0),find(timing_Lresponse > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
                                   Rr = intersect(find(timing_Rresponse < time0),find(timing_Rresponse > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
                                   
                                   if (length(La) == 1 && length(Lr) == 1) || (length(Ra) == 1 && length(Rr) == 1) % correct response
                                      onsets{4} = [onsets{4}; time0];
                                   end
                                   
                               end
                               
                               for i = 1:length(timing_feedback_2SP)
                                   time0 = timing_feedback_2SP(i);
                                   La = intersect(find(timing_Larrow < time0),find(timing_Larrow > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
                                   Ra = intersect(find(timing_Rarrow < time0),find(timing_Rarrow > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
                                   Lr = intersect(find(timing_Lresponse < time0),find(timing_Lresponse > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
                                   Rr = intersect(find(timing_Rresponse < time0),find(timing_Rresponse > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
                                   
                                   if (length(La) == 1 && length(Lr) == 1) || (length(Ra) == 1 && length(Rr) == 1) % correct response
                                      % now check if this is a win trial or a no win trial
                                      onsets{6} = [onsets{6}; time0];
                                   end
                                   
                               end
                               
                               for i = 1:length(timing_feedback_0SP)
                                   time0 = timing_feedback_0SP(i);
                                   La = intersect(find(timing_Larrow < time0),find(timing_Larrow > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
                                   Ra = intersect(find(timing_Rarrow < time0),find(timing_Rarrow > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
                                   Lr = intersect(find(timing_Lresponse < time0),find(timing_Lresponse > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
                                   Rr = intersect(find(timing_Rresponse < time0),find(timing_Rresponse > time0-2.0)); % actually -1.5s but due to small differences we can miss one. 
                                   if (length(La) == 1 && length(Lr) == 1) || (length(Ra) == 1 && length(Rr) == 1) % correct response
                                      % now check if this is a win trial or a no win trial
                                      index_trial = find(abs(timing_feedback-time0) == min(abs(timing_feedback-time0)));
                                      if order_trials(index_trial) == 2
                                         onsets{5} = [onsets{5}; time0];
                                      elseif order_trials(index_trial) == 4
                                         onsets{7} = [onsets{7}; time0];
                                      elseif order_trials(index_trial) == 5
                                         onsets{8} = [onsets{8}; time0];
                                      else
                                         fprintf('Something is wrong: inconsistencies found for feedback_C0 \n');
                                      end
                                   elseif (length(La) == 1 && length(Rr) == 1) || (length(Ra) == 1 && length(Lr) == 1)  % incorrect response
                                      % now check if this is a win trial or a no win trial
                                      index_trial = find(abs(timing_feedback-time0) == min(abs(timing_feedback-time0)));
                                      if order_trials(index_trial) == 2
                                         onsets{5} = [onsets{5}; time0];
                                      elseif order_trials(index_trial) == 4
                                         onsets{7} = [onsets{7}; time0];
                                      elseif order_trials(index_trial) == 5
                                         onsets{8} = [onsets{8}; time0];
                                      else
                                         fprintf('Something is wrong: inconsistencies found for feedback_C0 \n');
                                      end
                                      
                                   end
                                   
                               end
                           
                           % conversion to events.tsv table format - Lukas' code
                           
                           onset = [];
                           trial_type = categorical();
                            
                            for cond = 1:size(FID_trial_types,2)
                                
                                    trial_type_event = categorical();

                                    for ev = 1:size(onsets{cond},1)
                                        trial_type_event(ev,1) = FID_trial_types{cond};
                                    end

                                    onset = [onset; onsets{cond}];
                                    trial_type = [trial_type; trial_type_event];
                                    clear trial_type_event

                            end
                                
                           log=table(trial_type,onset);
                           log.duration = zeros(height(log),1);
                           log=sortrows(log,'onset','ascend');
                             
                           filename = fullfile(sessubjBIDSdir,[sourcesubjs{sub},['_ses-' num2str(ses)],'_task-',taskname2,runnames{run},'_events.tsv']);
                           writetable(log,filename,'Filetype','text','Delimiter','\t');
                           clear datafile index_* tmp* log onset* duration* trial_type filename
                           
                        end % if loop determining run number

                    end % if loop checking whether logfile exists

                end % for loop runs
                
                pheno_file_subj = [pheno_file_subj;pheno_file_ses];
                clear pheno_file_ses;
            
            end % for loop sessions
                        
        pheno_file = [pheno_file;pheno_file_subj];
        clear pheno_file_subj;

        end % for loop subjects
        
        pheno_filename = fullfile(pheno_dir,pheno_name);
        writetable(pheno_file,pheno_filename,'Filetype','text','Delimiter','\t');
    
    end % if loop checking sourcesubjs == BIDSsubjs
    
end % if loop checking writing option