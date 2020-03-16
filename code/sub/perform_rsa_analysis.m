function perform_rsa_analysis(p)
%
% - This code is written for performing multiple ROI representational similarity analysis.
% - A total of 383 (HCP360, 13 individual ROIs, and 10 subcortical ROIs) are used.
%
%
% called from emotion2020_analysis_BATCH
%
%% settinig
warning off

%% start analysis
for scoritr = 1:length(p.scoreTypes)
    scoreType = p.scoreTypes{scoritr};
    
    % RSA is performend only with four fature types
    switch scoreType
        case {'dimension','category','semantic','vision'};
        otherwise
            continue
    end
    
    for sbjitr = 1:length(p.sbjID)
        sbjID = p.sbjID{sbjitr};
        
        cntroi = 0;
        for roitr = randsample(1:length(p.roiDescrip),length(p.roiDescrip))%1:length(p.roiDescrip)
            cntroi = cntroi+1;
            roi = p.roiDescrip{roitr};
            tic
            suffix_full = [sbjID,'/',p.suffix_rsa,'/',roi];
            
            %% Save info. for final results
            saveFnameChk         = sprintf('%s/%s/%s/_log.txt',p.logdir,scoreType,suffix_full); % log files
            saveFname            = sprintf('%s/%s/%s/res.mat',p.savdir,scoreType,suffix_full); % res files
            setdir(fileparts(saveFnameChk));
            setdir(fileparts(saveFname));
            
            if p.checkModeRes
                chkfile = saveFname;
            else
                chkfile = saveFnameChk;
            end
            if exist(chkfile,'file')
                if exist(saveFname,'file')
                else
                    if p.del
                        fprintf('Delete%s\n',saveFnameChk)
                        delete(saveFnameChk)
                    end
                end
            elseif ~p.del
                fprintf('Start:%s\n',saveFnameChk)
                fprintf('CV %s [%s][%d/%d]:\n',suffix_full,scoreType,cntroi,length(p.roiDescrip))
                tmp = [];
                save(saveFnameChk,'tmp', '-ascii') % save log file
                
                % =======================
                % load data and get params for all subjects
                fprintf('Load data(%s)...\n',sbjID)
                dpath = sprintf('%s%s/rois/%s_%s.mat',p.fmridir,sbjID,sbjID,roi);
                load(dpath,'braindat','metainf');
                [nSample,nVox] = size(braindat);
                inds_all = 1:nSample;
                label_index = metainf.Label;
                inds_all = inds_all(~ismember(label_index,p.dupidx)); % ignore duplicate samples
                
                % load features
                load([p.featdir,scoreType,'.mat'],'L');
                L.feat = L.feat(label_index,:); % sort in presented order
                
                % z normalized
                braindat = zscore(braindat);
                L.feat = zscore(L.feat);
                
                clear rsa res
                brain_corr = getTri(fcorr(braindat(inds_all,:)'));
                feat_corr = getTri(fcorr(L.feat(inds_all,:)'));
                res.corr = corr(brain_corr,feat_corr);
                rsa.brain_corr = single(brain_corr);
                
                % display cs and accuracy
                fprintf('Representational similarity.\n')
                fprintf('rs = %.3f\n',nanmean(res.corr))
                tims
                
                fprintf('%s\n',saveFname)
                save(saveFname,'res','rsa','-v7.3')
                clear res rsa
            end
        end
    end
end

%%
end % end function
