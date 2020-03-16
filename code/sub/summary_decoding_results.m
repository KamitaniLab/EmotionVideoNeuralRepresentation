function decRes = summary_decoding_results(p)
%
% - This code is written for summarizing multiple ROI decoding analysis.
%
%
% called from emotion2020_analysis_BATCH
%
%% setting for ensemble decoders
thval = 0.095; % threshold for ROI selection
nAllrois = length(p.roiDescrip);

%% construct ensemble decoder

for scoritr = 1:length(p.scoreTypes)
    scoreType = p.scoreTypes{scoritr};
    
    % decoding is performend only with four fature types
    switch scoreType
        case {'dimension','category','semantic','vision'};
        otherwise
            continue
    end
    
    for sbjitr = 1:length(p.sbjID)
        sbjID = p.sbjID{sbjitr};
        
        % load metainf
        dpath = sprintf('%s%s/rois/%s_WholeBrain.mat',p.fmridir,sbjID,sbjID);
        load(dpath,'metainf');
        label_index = metainf.Label;
        
        suffix_summary = [sbjID,'/',p.suffix_dec,'/summary',];
        
        % save data inf.
        saveFname = sprintf('%s/%s/%s/res_mROIs_summary.mat',p.savdir,scoreType,suffix_summary);
        saveFname_feat = sprintf('%s/%s/%s/res_mROIs_feat.mat',p.savdir,scoreType,suffix_summary); % feat for UMAP analysis
        saveFnameChk = sprintf('%s/%s/%s/res_mROIs_summary_log.txt',p.logdir,scoreType,suffix_summary);
        setdir(fileparts(saveFname));
        setdir(fileparts(saveFnameChk));
        
        if exist(saveFnameChk,'file')
            fprintf('Load:%s\n',saveFname)
            load(saveFname,'ensDec','mRoiDec')
        else
            fprintf('Prepare ensemble:%s\n',saveFnameChk)
            tmp = [];
            save(saveFnameChk,'tmp','-ascii')
            
            % load features
            load([p.featdir,scoreType,'.mat'],'L');
            L.feat = L.feat(label_index,:); % sort in presented order
            nFeature = size(L.feat,2);
            L.feat = zscore(L.feat);
            
            
            % prepare for separate training and test samples
            nSample = size(L.feat,1);
            inds_all = 1:nSample;
            inds_all = inds_all(~ismember(label_index,p.dupidx));
            
            % initialization
            profile_cvacc_all = cell(nAllrois,1);
            profile_cvmxidx_all = cell(nAllrois,1);
            profile_acc_all = cell(nAllrois,1);
            pattern_acc_all = cell(nAllrois,1);
            iden_acc_all = cell(nAllrois,1);
            ensemble_select_idx = false(nAllrois,1); % whether within selected ROIs
            ensemble_accth_idx = false(nAllrois,nFeature); % whehter cv perf beyond threshold
            
            % load roi results and evalute performance while preparing ensemble
            clear ensDec
            ensDec.pred = cell(nFeature,1);
            ensDec.cnt = zeros(nFeature,1);
            ensDec.bestROIpred = cell(nFeature,1);
            ensDec.bestROIacc = zeros(nFeature,1);
            for roitr = 1:nAllrois
                roi = p.roiDescrip{roitr};
                
                suffix_full = [sbjID,'/',p.suffix_dec,'/',roi];
                dataFname = sprintf('%s/%s/%s/res.mat',p.savdir,scoreType,suffix_full);
                tmp = load(dataFname,'cv','r_best','iden');
                
                % evaluate performance by best model
                fprintf('[%03d]:Evaluate performance of %s:%s\n',roitr,scoreType,roi)
                acc_profile = tmp.r_best.all_best;
                acc_patttern = tmp.r_best.all_pattern_best;
                iden_acc = tmp.iden.cr*100;
                fprintf('mean profile acc.:%.3f\n',mean(acc_profile))
                fprintf('mean pattern acc.:%.3f\n',mean(acc_patttern))
                fprintf('mean identif acc.:%.1f%%\n',mean(iden_acc))
                
                profile_acc_all{roitr} = acc_profile;
                pattern_acc_all{roitr} = acc_patttern;
                iden_acc_all{roitr} = iden_acc;
                profile_cvacc_all{roitr} = tmp.cv.cvScoreAccEach_feat_bestVal;
                profile_cvmxidx_all{roitr} = tmp.cv.cvScoreAccEach_feat_bestIdx;
                
                % reserve prediction for ensemble decoder
                % select based roiset
                if sum(ismember(p.roiEnsemble,roi)) == 0
                    fprintf('Omitted from ensemble [not selected].\n')
                    continue
                else
                    ensemble_select_idx(roitr) = 1;
                    % select based on accuracy
                    for unitr = 1:nFeature
                        if tmp.cv.cvScoreAccEach_feat_bestVal(unitr) >= thval
                            ensemble_accth_idx(roitr,unitr) = 1;
                            if ensDec.cnt(unitr) == 0 % first roi
                                ensDec.pred{unitr} = tmp.r_best.pred_best(:,unitr);
                                ensDec.cnt(unitr) = 1;
                                ensDec.bestROIpred{unitr} = tmp.r_best.pred_best(:,unitr);
                                ensDec.bestROIacc(unitr) = profile_acc_all{roitr}(unitr);
                            else
                                ensDec.pred{unitr} = ensDec.pred{unitr}+tmp.r_best.pred_best(:,unitr);
                                ensDec.cnt(unitr) = ensDec.cnt(unitr)+1;
                            end
                        else
                            %fprintf('Omitted from ensemble [poor acc].\n')
                        end
                        % update best?
                        if profile_acc_all{roitr}(unitr) > ensDec.bestROIacc(unitr); % better for new.
                            ensDec.bestROIpred{unitr} = tmp.r_best.pred_best(:,unitr);
                            ensDec.bestROIacc(unitr) = profile_acc_all{roitr}(unitr);
                        end
                    end
                end
                tims
            end
            
            %% finalize ensemble
            ensDec.pred_final = zeros(length(inds_all),nFeature);
            for unitr = 1:nFeature
                if ensDec.cnt(unitr) == 0
                    fprintf('Max ROI was used as ensemble: %s:%s:%s\n',sbjID,scoreType,L.featname{unitr})
                    ensDec.pred_final(:,unitr) = ensDec.bestROIpred{unitr};
                else
                    ensDec.pred_final(:,unitr) = ensDec.pred{unitr}./ensDec.cnt(unitr);
                end
            end
            
            fprintf('Ensemble:\n')
            ensDec.profile_acc = diag(fcorr(ensDec.pred_final,L.feat(inds_all,:)));
            fprintf('Profile acc:%.3f\n',mean(ensDec.profile_acc))
            
            fprintf('Identification for integrated ROI\n')
            nanFeatIdx = isnan(ensDec.profile_acc);
            fc = fcorr(ensDec.pred_final(:,~nanFeatIdx)',L.feat(inds_all,~nanFeatIdx)');
            nanIdx = any(isnan(fc));
            ensDec.pred_true_corr = single(fc);
            [ensDec.idenacc_mean(~nanIdx,1),...
                ci,ensDec.idenacc_wins] = pwidentification(fc(~nanIdx,~nanIdx),1:(length(inds_all)-sum(nanIdx)));
            
            fprintf('Iden acc[%s][%d/%d]: cr = %.3f\n',scoreType,sum(~nanFeatIdx),length(nanFeatIdx),...
                nanmean(ensDec.idenacc_mean))
            
            clear mRoiDec
            mRoiDec.profile_acc_all = single(merge(profile_acc_all,2));
            mRoiDec.profile_cvacc_all = single(merge(profile_cvacc_all,2));
            mRoiDec.profile_cvmxidx_all = single(merge(profile_cvmxidx_all,2));
            mRoiDec.pattern_acc_all = single(merge(pattern_acc_all,2));
            mRoiDec.iden_acc_all = single(merge(iden_acc_all,2));
            mRoiDec.roi_all = p.roiDescrip;
            mRoiDec.roi_ensemble = p.roiEnsemble';
            mRoiDec.ensemble_select_idx = ensemble_select_idx;
            mRoiDec.ensemble_accth_idx = ensemble_accth_idx;
            
            fprintf('Save:%s\n',saveFname)
            save(saveFname,'ensDec','mRoiDec','-v7.3')
            
            feat_decoded = ensDec.pred_final;
            feat_true = L.feat(inds_all,:);
            save(saveFname_feat,'feat_true','feat_decoded','-v7.3')
            
        end
        
        decRes.(scoreType){sbjitr}.ensDec = ensDec;
        decRes.(scoreType){sbjitr}.mRoiDec = mRoiDec;
    end
end

%%
