function perform_emotion_identification_analysis(p,decRes)
%
% - This code is written for performing emotion identification analysis
%   based on region-wise decoding accuracy
%
%
% called from emotion2020_analysis_BATCH
%

%% settings
% score index
dimIdx = 1:14;
catIdx = 15:48;

hcpIdx = ~cellfun(@isempty,strfind(p.roiDescrip,'hcp180'));
subIdx = ismember(p.roiDescrip,{'thalamus','hyppocampus','hypothalamus','pallidum','brainstem','caudate','putamen','nuc_accum','amygdala','cerebellum'});
roiIdx = hcpIdx|subIdx;

nSbj = length(p.sbjID);


%% simiarlity matrix of various models

% identification parameters
nrepetitions = 100;
ncandidates = 2:length(catIdx);
ncandidatesdim = 2:length(dimIdx);

saveFname = sprintf('%s/emotion_identification/result_%d.mat',p.savdir,nrepetitions);
saveFnameChk = sprintf('%s/emotion_identification/result_%d_log.txt',p.savdir,nrepetitions);
setdir(fileparts(saveFname));
setdir(fileparts(saveFnameChk));

if ~exist(saveFnameChk,'file')
    tmp = [];
    save(saveFnameChk,'tmp','-ascii')
    
    DimAccROIs = cell(length(p.sbjID),1);
    CatAccROIs = cell(length(p.sbjID),1);
    for sbjitr = 1:length(p.sbjID)
        dimRes = decRes.dimension{sbjitr};
        catRes = decRes.category{sbjitr};
        
        DimAccROIs{sbjitr} = dimRes.mRoiDec.profile_acc_all(:,roiIdx)';
        CatAccROIs{sbjitr} = catRes.mRoiDec.profile_acc_all(:,roiIdx)';
        
    end
    
    % identification based on ROI-pattern similarity across individuals
    CatSimMat = cell(nSbj);
    CatIdenAcc = cell(nSbj);
    CatIdenAcc_mu = cell(nSbj);
    DimSimMat = cell(nSbj);
    DimIdenAcc = cell(nSbj);
    DimIdenAcc_mu = cell(nSbj);
    for sbjitr1 = 1:nSbj
        for sbjitr2 = 1:nSbj
            fprintf('%d:%d\n',sbjitr1,sbjitr2)
            if sbjitr1 == sbjitr2
                continue
            end
            CatSimMat{sbjitr1,sbjitr2} = fcorr(CatAccROIs{sbjitr1},CatAccROIs{sbjitr2});
            DimSimMat{sbjitr1,sbjitr2} = fcorr(DimAccROIs{sbjitr1},DimAccROIs{sbjitr2});
            
            CatIdenAcc{sbjitr1,sbjitr2} = mcidentification(CatSimMat{sbjitr1,sbjitr2},1:length(catIdx),ncandidates,nrepetitions);
            DimIdenAcc{sbjitr1,sbjitr2} = mcidentification(DimSimMat{sbjitr1,sbjitr2},1:length(dimIdx),ncandidatesdim,nrepetitions);
            CatIdenAcc_mu{sbjitr1,sbjitr2} = cellfun(@mean,CatIdenAcc{sbjitr1,sbjitr2});
            DimIdenAcc_mu{sbjitr1,sbjitr2} = cellfun(@mean,DimIdenAcc{sbjitr1,sbjitr2});
        end
    end
    
    CatIdenAcc_all = merge([getTri(CatIdenAcc);getTri(CatIdenAcc')],2);
    CatSimMat_mu = cellmean([getTri(CatSimMat);getTri(CatSimMat')]);
    CatSimMat_all = ([getTri(CatSimMat);getTri(CatSimMat')]);
    DimIdenAcc_all = merge([getTri(DimIdenAcc);getTri(DimIdenAcc')],2);
    DimSimMat_mu = cellmean([getTri(DimSimMat);getTri(DimSimMat')]);
    DimSimMat_all = ([getTri(DimSimMat);getTri(DimSimMat')]);
    
    [CatIdenAcc_ci_all, CatIdenAcc_mu_all] = ciestim3(merge([getTri(CatIdenAcc_mu);getTri(CatIdenAcc_mu')],2)*100,2);
    [DimIdenAcc_ci_all, DimIdenAcc_mu_all] = ciestim3(merge([getTri(DimIdenAcc_mu);getTri(DimIdenAcc_mu')],2)*100,2);
    
    fprintf('Save:%s\n',saveFname)
    save(saveFname,'CatIdenAcc_mu_all','CatIdenAcc_ci_all','DimIdenAcc_mu_all','DimIdenAcc_ci_all',...
        'CatIdenAcc_all','DimIdenAcc_all','CatIdenAcc_mu','DimIdenAcc_mu','CatSimMat_all','DimSimMat_all',...
        '-v7.3');
end

%%
