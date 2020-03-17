% function emotion2020_analysis_BATCH
%
% This code include
%   - Data preparations for delineating indivual ROIs, including WholeBrain, HCP360 ROIs, subcortical regions
%   - (regularized) linear regression analyses (encoding/decoding) between MRI data 
%      and labels (category, dimension, vision, semantic) associated with 2196 (2181 unique) emotion evocative movie clips.
%   - Representational similarity analysis
%   - K-means clustering using emotion-related brain activity (encoding results are necessary)
%
% Preparations:
%   - Set preprocessed fmri data in root/data/fmri/SubjectX/preprocessed/ 
%   - Set roiInf.mat file in root/code/data/fmri/misc/
%   - Set feature data in root/data/features/ 
%   - Set principal gradient data in root/data/fmri/SubjectX/pringrad/ 
%   - Set BrainDecodeeerToolbox2 in root/code/libraries/ 
%
% Main parts:
%   - To go through all analyses and get all result figures, run this and python scripts as below.
%     0. run this script with setting 1 for roiDataPreparation variable (multiple cpu can work in parallel)
%     1. run this script with setting 1 for performDecAnalyses/performEncAnalyses/performRSAnalyses variables (multiple cpu can work in parallel)
%     2. run this script with setting 1 for summaryDecAnalyses/summaryEncAnalyses/summaryRSAnalyses variables
%     3. run this script with setting 1 for performeAdditionalAnalysis variable
%     4. run python scripts to perform UMAP analyses.
%     5. run this script with setting 1 for showDecResults/showEncResults/showRSAResults variable
%
% Note:
%   - Decoding, encoding, and representational similarity analyses can be performed independently.
%   - Two umap analysis implemented in python scripts requires results of decoding and encoding analyses.
%   - The analysis part (2) will take about 1 day using 100 cpu to complete all the computations.
%
%
% written by Tomoyasu Horikawa horikawa-t@atr.jp 2020/03/11
%

%% initialize
clear all, close all

% general info.
p.rootPath = './'; % set the current directory
p.fileName = 'emotion2020_analysis_BATCH';

%% analyses setting [CHANGE this section for switching]
roiDataPreparation = 0; % set 1 to prepare individual ROI data

performDecAnalyses = 0; % set 1 to perform decoding analyses
performEncAnalyses = 0; % set 1 to perform encoding analyses
performRSAnalyses  = 0; % set 1 to perform representational similarity analyses

summaryDecAnalyses = 0; % set 1 to summarize decoding analyses
summaryEncAnalyses = 0; % set 1 to summarize encoding analyses
summaryRSAnalyses  = 0; % set 1 to summarize representational similarity analyses

performeAdditionalAnalysis = 0; % set 1 to perform additional analyses (emotion identification, clustering)

showDecResults = 0; % set 1 to show results of decoding analyses
showEncResults = 0; % set 1 to show results of encoding analyses
showRSAResults = 0; % set 1 to show results of representational similarity analyses

% misc setting
p.del = 0; % set 1 if you want to clean log files (delete log files without res files)
p.checkChkfile = 0; % if 0 skip file check (fast) .
p.checkModeRes = 0; % if 1 check by result file, if 0 check by log file.
p.integrateRes = 1; % if 1 integrate results across nparse

% add path to code
addpath(genpath([p.rootPath,'code/']));
pathset([p.rootPath,'code/libraries/BrainDecoderToolbox2/fig/'],0); % to avoid duplicates

%% set path information
% set directories
p.savdir = setdir(sprintf('%s/res/%s/',p.rootPath,p.fileName)); % analysis results
p.logdir = setdir(sprintf('%s/log/%s/',p.rootPath,p.fileName)); % log files
p.figdir = setdir(sprintf('%s/fig/',p.rootPath)); % result figures

p.fmridir = setdir(sprintf('%s/data/fmri/',p.rootPath));
p.miscdir = setdir(sprintf('%s/data/fmri/misc/',p.rootPath));
p.featdir = setdir(sprintf('%s/data/feature/',p.rootPath));
p.roiInfPath = [p.miscdir,'roiInf.mat'];

%% data information
% subject info
p.sbjID = {'Subject1','Subject2','Subject3','Subject4','Subject5'};

% duplicate index (label index of duplicated videos). These wll be removed.
p.dupidx = [1,4:8,11,859,866,1673,2157,2187,2188,2194,2195];

%% label group types
p.scoreTypes = {'category','dimension','dim_cat','semantic','vision','dim28binary','dim28continuous','categcontinuous'};

%% analysis parameters
% cross-validation parameters
p.run2FoldAssignIdx = [1,11,21,31,41,51;
    10,20,30,40,50,61]; % 1st/2nd row = start/end run indices for each CV fold.
p.nCVfolds = size(p.run2FoldAssignIdx,2);

% reguralization paramters settings
p.nparamLogSearch = 20; %  You can reduce computation time by setting smaller # of exploration parameters.
p.lambda = logspace(1,4,p.nparamLogSearch);

% decoding analysis parameters
p.nSelectVoxels = 500; % # of selected voxels (n = 500 in our study). You can reduce computation time by setting smaller # of voxels (e.g., 50) here.
p.nparse = 100; % number of groups for whole brain encoding analyses

%% set ROIs for ROIbased analysis (decoding, rsa)
load(p.roiInfPath,'roiInf');
p.roiDescrip = {roiInf.hcp360.roiDescrip{:},roiInf.individual.roiDescrip{:},roiInf.subcortex.roiDescrip{:}};
p.roiEnsemble = {roiInf.hcp360.roiDescrip{:},roiInf.subcortex.roiDescrip{:}}; % roi sets for ensemble decoders

%% Preparation section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% data prepartion (must be done before starting analyses)
if roiDataPreparation
    prepare_roidata(p.rootPath)
end

%% Analyses section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% perform decoding ROI based analyses
p.suffix_dec = ['nCV',num2str(p.nCVfolds),'_mROIDec/'];
if performDecAnalyses
    perform_decoding_analysis(p)
end

%% perform encoding whole brain analyses
p.suffix_enc = ['nCV',num2str(p.nCVfolds),'_WholeBrainEnc/'];
if performEncAnalyses
    perform_encoding_analysis(p)
end

%% perform representational similarity analyses
p.suffix_rsa = ['_mROIRSA/'];
if performRSAnalyses
    perform_rsa_analysis(p)
end

%% Load/summary result section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if summaryDecAnalyses
    decRes = summary_decoding_results(p);
end
if summaryEncAnalyses
    encRes = summary_encoding_results(p);
end
if summaryRSAnalyses
    rsaRes = summary_rsa_results(p);
end

%% additional analysis section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if performeAdditionalAnalysis
    perform_emotion_identification_analysis(p,decRes)
    perform_clustering_analysis(p,encRes)
end

%% Draw result section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if showDecResults
    draw_figure2AB(p,decRes)
    draw_figure2DE(p,decRes)
    draw_figureS1AB(p,decRes)
    draw_figure3ABC_S1C(p,decRes)
    draw_figure3D(p,decRes) % umap analysis need to be done beforehand
end
if showEncResults
    draw_figure4BC_S2B(p,encRes)
    draw_figureS3AB(p,encRes)
    draw_figure4D_S4D(p,encRes)
    draw_figure4E_5BCD(p,encRes)
    draw_figure4F(p,encRes)
    draw_figure5GH(p,encRes)
    draw_figure6AB(p) % umap analysis need to be done beforehand
    draw_figure6CDE_S6(p)
end
if showDecResults && showRSAResults
    draw_figureS5(p,decRes,rsaRes)
end


%%
