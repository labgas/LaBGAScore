%% LaBGAScore_prep_graphvar_input_from_conn.m
%
%
% *USAGE*
%
% This script creates the necessary input files for graph theoretical
% analysis with the GraphVar toolbox from ROI-to-ROI connectivity analysis
% performed with the CONN toolbox
%
% NOTE: this script is not based on LaBGAS standard data structure for fmri
%       projects, but can easily be adapted/simplified in the future
%
% NOTE: make a new GraphVar workspace using the GraphVar GUI prior to
%       running this script, and enter its name in the option below!
%
% Launch GraphVar from your Matlab terminal by running start_GraphVar
%
% For more info about GraphVar, see the following resources
%
% # GraphVar website: http://rfmri.org/GraphVar
%
% # GraphVar YouTube channel: https://www.youtube.com/@graphvar9022
%
%
% The script also includes an option to generate input files for the
% NBS-predict and BrainNetClass toolboxes
%
% For more info about NBS-predict, see the following resources
%
% # NBS-predict Github page: https://github.com/eminSerin/NBS-Predict
%
% # NBS-predict YouTube tutorial: https://www.youtube.com/watch?v=PiENWTFCLUo
%
% For more info about BrainNetClass, see the following resources
%
% # BrainNetClass Github page: https://github.com/zzstefan/BrainNetClass
%
% # BrainNetClass paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7294070/
%
%
% *OPTIONS*
%
% * proj_name_conn          name of your CONN project
%
% * proj_name_graphvar      name of your GraphVar workspace
%                           NOTE: can be the same or different from CONN
%                           project, so you can make different GraphVar
%                           workspaces per CONN project (this is better
%                           than saving all interim results in one
%                           workspace)
%
% * input_data              choose input data type
%
%   1. 'corr_matrix'
%      Uses ROI-to-ROI correlation matrix calculated by CONN as input to 
%      GraphVar, and calculate parametric p-values from r-values and df
%      THIS OPTION WILL ALSO PRODUCE CORRELATION MATRICES FOR NBS-PREDICT
%
%   2. 'denoised_timeseries'
%      Uses denoised timeseries in each ROI as written by CONN as input for
%      ROI-to-ROI correlation matrix calculation by GraphVar, allowing for
%      more advanced options such as non-parametric p-value calculation, 
%      significance-based network thresholding, partial correlations, 
%      dynamic connectivity analysis (sliding windows), and fancy 
%      calculations as implemented in the BrainNetClass option in GraphVar
%      THIS OPTION WILL PRODUCE DENOISED TIMESERIES FILES FOR BRANNETCLASS
%
% * proj_name_nbspredict    (optional)
%                           if you want to generate input files for
%                           NBS-predict in addition to GraphVar input
%                           files, enter the project name in this option
%                           otherwise, comment out
%
% * proj_name_brainnetclass (optional)
%                           if you want to generate input files for
%                           brainnetclass in addition to GraphVar input
%                           files, enter the project name in this option
%                           otherwise, comment out
%
% 
% -------------------------------------------------------------------------
%
% author: Lukas Van Oudenhove
%
% date:   KU Leuven, July, 2023
%
% -------------------------------------------------------------------------
%
% LaBGAScore_prep_graphvar_input_from_conn.m         v1.2
%
% last modified: 2023/08/09
%
%
%% SET OPTIONS & PATHS
% -------------------------------------------------------------------------

% LaBGAScore_prep_s0_define_directories; % should work with standard LaBGAS
% datasets

proj_name_conn = 'proj-FD-PPI';
proj_name_graphvar = 'proj-FD-PPI_timeseries_corr';
phenofilename2load = 'phenotype.tsv';
phenofilename2write = 'Variables.xlsx';
input_data = 'denoised_timeseries';

projrootdir = 'C:\Users\lukas\Dropbox (Dartmouth College)\proj-FD-PPI';
conn_dirs = struct();
conn_dirs.connrootdir = fullfile(projrootdir,'analysis',['conn_' proj_name_conn]);
conn_dirs.conndatadir = fullfile(conn_dirs.connrootdir,'data');
conn_dirs.connresultsdir = fullfile(conn_dirs.connrootdir,'results');
conn_dirs.connpreprocdir = fullfile(conn_dirs.connresultsdir,'preprocessing');
conn_dirs.connfirstleveldir = fullfile(conn_dirs.connresultsdir,'firstlevel','RRC_01');

projBIDSdir = fullfile(projrootdir,'BIDS');
projphenodir = fullfile(projBIDSdir,'phenotype');

projderivdir = fullfile(projrootdir,'derivatives','fmriprep');
derivlist = dir(fullfile(projderivdir,'sub-*'));
derivlist = derivlist([derivlist(:).isdir]);
derivsubjs = cellstr(char(derivlist.name));
nr_derivsubjs = size(derivsubjs,1);

graphvarrootdir = 'C:\Users\lukas\Documents\MATLAB\GraphVar_2.03a';
graphvar_dirs = struct();
graphvar_dirs.graphvarworkspacesdir = fullfile(graphvarrootdir,'workspaces');
graphvar_dirs.graphvarprojdir = fullfile(graphvar_dirs.graphvarworkspacesdir, proj_name_graphvar);
graphvar_dirs.graphvardatadir = fullfile(graphvar_dirs.graphvarprojdir,'data');
graphvar_dirs.graphvarresultsdir = fullfile(graphvar_dirs.graphvarprojdir,'results');

conn_condition_name = 'Condition001';

% Optional for NBS-predict

proj_name_nbspredict = 'proj-FD-PPI';

    if exist('proj_name_nbspredict','var')
        nbspredict_dirs = struct();
        nbspredict_dirs.nbspredictrootdir = fullfile(projrootdir,'analysis',['NBSPredict_' proj_name_nbspredict]);
        nbspredict_dirs.nbspredictcorrmatdir = fullfile(nbspredict_dirs.nbspredictrootdir,'CorrMatrix');
            if ~exist(nbspredict_dirs.nbspredictcorrmatdir,'dir')
                mkdir(nbspredict_dirs.nbspredictcorrmatdir);
            end
    end
    
% Optional for BrainNetClass

proj_name_brainnetclass= 'proj-FD-PPI';

    if exist('proj_name_brainnetclass','var')
        brainnetclass_dirs = struct();
        brainnetclass_dirs.brainnetclassrootdir = fullfile(projrootdir,'analysis',['BrainNetClass_' proj_name_brainnetclass]);
        brainnetclass_dirs.brainnetclasstimeseriesdir = fullfile(brainnetclass_dirs.brainnetclassrootdir,'TimeSeries');
            if ~exist(brainnetclass_dirs.brainnetclasstimeseriesdir,'dir')
                mkdir(brainnetclass_dirs.brainnetclasstimeseriesdir);
            end
        brainnetclass_dirs.brainnetclassresultsdir = fullfile(brainnetclass_dirs.brainnetclassrootdir,'Results');
            if ~exist(brainnetclass_dirs.brainnetclassresultsdir,'dir')
                mkdir(brainnetclass_dirs.brainnetclassresultsdir);
            end
    end
    

%% CREATE GRAPHVAR BRAIN REGIONS FILE
% -------------------------------------------------------------------------

% NOTE: this code is atlas-specific, this is for JHU atlas

brain_regions_table = table();
connsummary = load(fullfile(conn_dirs.connfirstleveldir,['resultsROI_' conn_condition_name '.mat']));
brain_regions_table.Var1 = ones(size(connsummary.names,2),1);

    for region = 1:height(brain_regions_table)
        name = strsplit(connsummary.names{region},' (');
        short_name = strsplit(name{1},'.');
        short_name = short_name{2};
        brain_regions_table.Var2(region) = {short_name};
        if size(name,2) == 2
            long_name = name{2}(1:end-3);
        else
            long_name = [name{2} ' (' name{3}(1:end-3)];
        end
        brain_regions_table.Var3(region) = {long_name};
        clear name short_name long_name
        xyz = connsummary.xyz{region};
        brain_regions_table.Var4(region) = round(xyz(1,1));
        brain_regions_table.Var5(region) = round(xyz(1,2));
        brain_regions_table.Var6(region) = round(xyz(1,3));
        clear xyz
    end

writetable(brain_regions_table,fullfile(graphvar_dirs.graphvarprojdir,'BrainRegions.xlsx'),'WriteVariableNames',false);

    if exist('proj_name_nbspredict','var')
        brain_regions_table_nbspredict = brain_regions_table;
        brain_regions_table_nbspredict.Var7 = brain_regions_table_nbspredict.Var2;
        brain_regions_table_nbspredict = brain_regions_table_nbspredict(:,4:end);
        writetable(brain_regions_table_nbspredict,fullfile(nbspredict_dirs.nbspredictrootdir,'node_description.csv'),'WriteVariableNames',false);
    end


%% LOAD DATA FROM CONN, TRANSFORM, CALCULATE P-VALUES, AND SAVE
% -------------------------------------------------------------------------

switch input_data
    
    case 'corr_matrix'
        
        connmatlist = dir(fullfile(conn_dirs.connfirstleveldir,'resultsROI_Subject*'));
        connpreproclist = dir(fullfile(conn_dirs.connpreprocdir,'DATA_*.mat'));
        nr_connsubjs = size(connmatlist,1);
        
        if nr_connsubjs ~= nr_derivsubjs
            error('\nNumber of subjects in derivatives directory %s does not match number of subjects in conn directory %s, please check before proceeding/n', projderivdir, conn_dirs.connpreprocdir);
        end
        
        corr_matrices = cell(nr_connsubjs,1);
        pval_matrices = cell(nr_connsubjs,1);
        tval_matrices = cell(nr_connsubjs,1);
        connsubjs = cell(nr_connsubjs,1);
        
        for n = 1:nr_connsubjs
            
            connmat = load(fullfile(conn_dirs.connfirstleveldir,connmatlist(n).name));
            conn_subj = connmatlist(n).name(12:21);
            connsubjs{n} = conn_subj;
            
            load(fullfile(conn_dirs.connpreprocdir,connpreproclist(n).name));
            
            N = V.size.Nt;
            CorrMatrix = tanh(connmat.Z); % inverse Fisher r-to-z transform
            CorrMatrix(isnan(CorrMatrix))=1; % set diagonal elements to 1 rather than NaN
            
            TValMatrix = CorrMatrix.*sqrt((N-2)./(1-CorrMatrix.^2)); % convert r to t
            PValMatrix = 2.*(1-tcdf(abs(TValMatrix),(N-2))); % calculate p from t and df
            
            corr_matrices{n} = CorrMatrix;
            tval_matrices{n} = TValMatrix;
            pval_matrices{n} = PValMatrix;
            
            matfilename = fullfile(graphvar_dirs.graphvardatadir,'CorrMatrix',['CorrMatrix_' connsubjs{n} '.mat']);
            save(matfilename,'CorrMatrix','TValMatrix','PValMatrix','-v7.3');
            
            if exist('proj_name_nbspredict','var')
               matfilename_nbspredict = fullfile(nbspredictrootdir,'CorrMatrix',['CorrMatrix_' connsubjs{n} '.mat']);
               save(matfilename_nbspredict,'CorrMatrix','-v7.3'); % NBS-Predict cannot handle matfiles with multiple arrays
            end
            
            clear connmat conn_subj V N CorrMatrix TValMatrix PValMatrix matfilename;
            
        end
                
    case 'denoised_timeseries'
        
        connmatlist = dir(fullfile(conn_dirs.connpreprocdir,'ROI_Subject*_Condition001.mat'));
        nr_connsubjs = size(connmatlist,1);
        
        if nr_connsubjs ~= nr_derivsubjs
            error('\nNumber of subjects in derivatives directory %s does not match number of subjects in conn directory %s, please check before proceeding/n', projderivdir, conn_dirs.connpreprocdir);
        end
        
        roi_signals = cell(nr_connsubjs,1);
        connsubjs = cell(nr_connsubjs,1);
        
        for n = 1:nr_connsubjs
            
            load(fullfile(conn_dirs.connpreprocdir,connmatlist(n).name),'data');
            load(fullfile(conn_dirs.connpreprocdir,connmatlist(n).name),'names');
            conn_subj = connmatlist(n).name(5:14);
            connsubjs{n} = conn_subj;
            
            idx_names = contains(names,brain_regions_table.Var3);
            roi_signal = cell2mat(data(idx_names));
            
            roi_signals{n} = roi_signal;

%             matfilename = fullfile(graphvar_dirs.graphvardatadir,'Signals',['ROISignals_' connsubjs{n} '.mat']);
%             save(matfilename,'roi_signal','-v7.3');
            
            if exist('proj_name_brainnetclass','var')
               txtfilename_brainnetclass = fullfile(brainnetclass_dirs.brainnetclasstimeseriesdir,['ROISignals_' connsubjs{n} '.txt']);
               writematrix(roi_signal,txtfilename_brainnetclass);
            end
            
            clear conn_subj idx_names roi_signal matfilename;
            
        end
               
end


%% CREATE/LOAD/EDIT PHENOTYPE FILE AND WRITE GRAPHVAR VARIABLES FILE
% -------------------------------------------------------------------------

% create from scratch

pheno_table = table();
pheno_table.Subj_ID = connsubjs;
pheno_table.Subj_ID_derivs = derivsubjs;
for sub = 1:height(pheno_table)
    if contains(pheno_table.Subj_ID_derivs{sub},'patient')
        pheno_table.Group(sub,:) = 'FD';
    else
        pheno_table.Group(sub,:) = 'HC';
    end
end

writetable(pheno_table,fullfile(graphvar_dirs.graphvarprojdir,phenofilename2write));

    if exist('proj_name_nbspredict','var')
        design_matrix = zeros(height(pheno_table),2);
         for sub = 1:height(pheno_table)
             if strcmp(pheno_table.Group(sub,:),'HC')
                 design_matrix(sub,1) = 0;
                 design_matrix(sub,2) = 1;
             elseif strcmp(pheno_table.Group(sub,:),'FD')
                 design_matrix(sub,1) = 1;
                 design_matrix(sub,2) = 0;
             end
         end
        matfilename_nbspredict_design = fullfile(nbspredict_dirs.nbspredictrootdir,'design_matrix.mat');
        save(matfilename_nbspredict_design,'design_matrix','-v7.3');
    end
    
    if exist('proj_name_brainnetclass','var')
        labels = zeros(height(pheno_table),1);
         for sub = 1:height(pheno_table)
             if strcmp(pheno_table.Group(sub,:),'HC')
                 labels(sub,1) = -1;
             elseif strcmp(pheno_table.Group(sub,:),'FD')
                 labels(sub,1) = 1;
             end
         end
        txtfilename_brainnetclass_label = fullfile(brainnetclass_dirs.brainnetclassrootdir,'labels.txt');
        writematrix(labels,txtfilename_brainnetclass_label);
    end

% load and add CONN Subj_ID

% phenotable = readtable(fullfile(projphenodir,phenofilename2load),'Delimiter','t');
% pheno_table.Subj_ID = connsubjs;
% writetable(pheno_table,fullfile(graphvar_dirs.graphvarprojdir,phenofilename2write));
