function prepare_roidata(rootPath)
%  
% This script was written for preparing preprocessed MRI data for cowen &
% keltner movie clip presentation data, which were pre-processed by the fmriprep. 
% The data was devided into individual ROIs.
% 
% 
% call prepare_roidata_individual.m
% called from emotion2020_analysis_BATCH.m
% 
% Required libraries:
% BrainDecoderToolbox2
% 
% 
% Note: this script tentatively requires ~80GB memory for data preprocessing
% 
% written by Tomoyasu Horikawa horikawa-t@atr.jp 2020/03/11
% 

%% initialize
% general info.
if ~exist('rootPath','var')
    rootPath = '/home/nu/horikawa-t/toolbox/public/emotion2020/'; 
end
fileName = 'prepare_roidata';

% add path to code
addpath(genpath([rootPath,'code/']));

% set paths to directories/data
p.savdir = setdir(sprintf('%s/res/%s/',rootPath,fileName));
p.logdir = setdir(sprintf('%s/log/%s/',rootPath,fileName));
p.fmridir = setdir(sprintf('%s/data/fmri/',rootPath));
p.miscdir = setdir(sprintf('%s/data/fmri/misc/',rootPath));
p.roiInfPath = [p.miscdir,'roiInf.mat'];

%% subject info. {sbjID, nSession}
dataList = {
    {'Subject1',5}
    {'Subject2',7}
    {'Subject3',6}
    {'Subject4',5}
    {'Subject5',5}
    };
nSbj = length(dataList);

%% Set ROI information
p.roiDescrip = {};
p.roiName = {};
p.exRoiName = {};

% load prepared roiInf
load(p.roiInfPath,'roiInf');

% add individual ROIs for hcp360, 13 individual, and 10 subcortical ROIs
p.roiDescrip = {p.roiDescrip{:},roiInf.hcp360.roiDescrip{:},roiInf.individual.roiDescrip{:},roiInf.subcortex.roiDescrip{:}}';
p.roiName = {p.roiName{:},roiInf.hcp360.roiName{:},roiInf.individual.roiName{:},roiInf.subcortex.roiName{:}}';
p.exRoiName = {p.exRoiName{:},roiInf.hcp360.exRoiName{:},roiInf.individual.exRoiName{:},roiInf.subcortex.exRoiName{:}}';

% add whole brain ROIs
p.roiDescrip = {'WholeBrain',p.roiDescrip{:}}'; % cortex (HCP180) + specified subcortex
p.roiName = {[roiInf.hcp360.roiName{:},roiInf.subcortex.roiName{:}],p.roiName{:}}';
p.exRoiName = {{},p.exRoiName{:}}';

% add cortical ROIs
p.roiDescrip = {'Cortex',p.roiDescrip{:}}';
p.roiName = {[roiInf.hcp360.roiName{:}],p.roiName{:}}';
p.exRoiName = {{},p.exRoiName{:}}';

%% Preprocess MRI data
% subject loop
for datitr = 1:nSbj;
    % initialize parameter variable
    clear param
    param = p;
    
    % set path for individual subject
    param.sbj = dataList{datitr}{1};
    param.nSession = dataList{datitr}{2};
    param.fmridir_org = setdir([p.fmridir,param.sbj,'/preprocessed/']);
    param.fmridir_rois = setdir([p.fmridir,param.sbj,'/rois/']);
    param.logdir_sbj = setdir([p.logdir,param.sbj,'/']);
    
    % create data
    saveFname = [param.savdir,'dataPath_',param.sbj,'.mat'];
    saveFnameChk = [param.logdir_sbj,'dataPath_',param.sbj,'_log.txt'];
    if ~exist(saveFnameChk,'file')
        fprintf('Start: %s\n',saveFnameChk)
        tmp = [];
        save(saveFnameChk,'tmp','-ascii');
        param = prepare_roidata_individual(param);
        save(saveFname,'param','-v7.3')
    else
        fprintf('Data path file for %s: %s\n',param.sbj,saveFname)
    end
end

%%
