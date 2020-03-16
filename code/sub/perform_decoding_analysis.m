function perform_decoding_analysis(p)
%
% - This code is written for performing multiple ROI decoding analysis.
% - A total of 383 (HCP360, 13 individual ROIs, and 10 subcortical ROIs) are used.
%
%
% called from emotion2020_analysis_BATCH
%
%% initialization
warning off

%% start analysis
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
        
        cntroi = 0;
        for roitr = randsample(1:length(p.roiDescrip),length(p.roiDescrip))%1:length(p.roiDescrip)
            cntroi = cntroi+1;
            roi = p.roiDescrip{roitr};
            tic
            suffix_full = [sbjID,'/',p.suffix_dec,'/',roi];
            
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
                nFeature = size(L.feat,2);
                
                % z normalized
                braindat = zscore(braindat);
                L.feat = zscore(L.feat);
                
                % set containers for a voxel group
                clear res stats pred
                for pitr = 1:p.nparamLogSearch
                    res.scoreacc{pitr} = zeros(nFeature,1);
                    pred.prediction{pitr} = zeros(nFeature,1);
                    pred.prediction_nest{pitr} = zeros(nFeature,p.nCVfolds);
                    for ixx = 1:p.nCVfolds
                        res.scoreacc_nest{pitr,ixx} = zeros(nFeature,1);
                    end
                    res_tmp.tepred_score{pitr} = cell(1,p.nCVfolds);
                end
                
                %% analysis loop
                for ixx = 1:p.nCVfolds
                    % get sample index
                    inds_te = find(ismember(metainf.Run,p.run2FoldAssignIdx(1,ixx):p.run2FoldAssignIdx(2,ixx))); % get test index
                    label_index = metainf.Label(inds_te);
                    inds_te = inds_te(~ismember(label_index,p.dupidx));
                    inds_tr = inds_all(~ismember(inds_all,inds_te));
                    
                    clear D_tr D_te
                    % separate training and test data
                    D_tr.brain = braindat(inds_tr,:);
                    D_te.brain = braindat(inds_te,:);
                    D_tr.score = L.feat(inds_tr,:);
                    D_te.score = L.feat(inds_te,:);
                    
                    % nested CV ===============================
                    % nested CV was done within the Training data
                    clear res_tmp_nest
                    trainRunIdx = ~ismember(1:p.nCVfolds,ixx);
                    nestcvIdx = [p.run2FoldAssignIdx(1,trainRunIdx);p.run2FoldAssignIdx(2,trainRunIdx)];
                    p.nCVfolds_nest = size(nestcvIdx,2);
                    Run_train = metainf.Run(inds_tr);
                    inds_all_nest = 1:length(inds_tr);
                    
                    % container for each voxel.
                    for pitr = 1:p.nparamLogSearch
                        res_tmp_nest.tepred_score{pitr,ixx} = cell(1,p.nCVfolds_nest);
                    end
                    
                    for ixx2 = 1:p.nCVfolds_nest
                        clear D_tr_nest D_te_nest
                        % separate training and test data
                        inds_te_nest = find(ismember(Run_train,nestcvIdx(1,ixx2):nestcvIdx(2,ixx2))); % get test index
                        inds_tr_nest = find(~ismember(inds_all_nest,inds_te_nest)); % get traiing index
                        D_tr_nest.brain = D_tr.brain(inds_tr_nest,:);
                        D_te_nest.brain = D_tr.brain(inds_te_nest,:);
                        D_tr_nest.score = D_tr.score(inds_tr_nest,:);
                        D_te_nest.score = D_tr.score(inds_te_nest,:);
                        
                        % brain normalization
                        parm.data_norm = 1;
                        [train_brain,nparm] = normalize_data(D_tr_nest.brain',parm.data_norm);
                        parm.xmean = nparm.xmean;
                        parm.xnorm = nparm.xnorm;
                        test_brain = normalize_data(D_te_nest.brain',parm.data_norm,parm);
                        
                        % score normalization
                        [train_score,nparm] = normalize_data(double(D_tr_nest.score)',parm.data_norm);
                        
                        % perform l2 analysis
                        % voxel selection
                        cor = fcorr(train_brain',train_score');
                        nanIdx = any(isnan(cor),2);
                        cor(isnan(cor(:))) = 0; % set 0 if nan
                        [sorted,ord] = sort(abs(cor),'descend');
                        useVoxCand = sort(ord(1:min([p.nSelectVoxels,nVox]),:),1);
                        [uniUseVoxCand] = unique(useVoxCand(:));
                        uniUseVoxCand(ismember(uniUseVoxCand,find(nanIdx))) = [];
                        train_brain_use = train_brain(uniUseVoxCand,:)';
                        test_brain_use = test_brain(uniUseVoxCand,:);
                        
                        % add bias term
                        Xtr_1 = addBias(train_brain_use);
                        % calcurate weights
                        X1X1 = Xtr_1'*Xtr_1;
                        X1T = Xtr_1'*train_score';
                        for i = 1:nFeature
                            useVoxIdx = find(ismember(uniUseVoxCand,useVoxCand(:,i)));
                            for pitr = 1:p.nparamLogSearch
                                Lam = eye(length(useVoxIdx)+1)*p.lambda(pitr);
                                XXL = X1X1([useVoxIdx;end],[useVoxIdx;end])+Lam;
                                W = XXL\X1T([useVoxIdx;end],i);
                                pred_tmp = addBias(test_brain_use(useVoxIdx,:)')*W;
                                res_tmp_nest.tepred_score{pitr,ixx}{ixx2}(:,i) = pred_tmp.*nparm.xnorm(i) + nparm.xmean(i);
                            end
                        end
                    end
                    
                    % decoding summary in nestCV
                    for pitr = 1:p.nparamLogSearch
                        tepred_score = merge(res_tmp_nest.tepred_score{pitr,ixx},1); % merge prediction
                        scoreacc = fcorrdiag(D_tr.score,tepred_score); % evaluate correlations between true and pred
                        res.scoreacc_nest{pitr,ixx} = single(scoreacc); % reserve prediction accuracy
                    end
                    clear D_tr_nest D_te_nest parm parm_all_nest train_score test_score train_brain test_brain res_tmp_nest
                    % end nested CV analyses
                    % ========================
                    
                    
                    % res-start CV analyses (not nestedCV)========================
                    parm.data_norm = 1;
                    [train_brain,nparm] = normalize_data(D_tr.brain',parm.data_norm);
                    parm.xmean = nparm.xmean;
                    parm.xnorm = nparm.xnorm;
                    test_brain = normalize_data(D_te.brain',parm.data_norm,parm);
                    
                    % score normalization
                    [train_score,nparm] = normalize_data(double(D_tr.score)',parm.data_norm);
                    
                    % perform l2 analysis
                    % voxel selection
                    cor = fcorr(train_brain',train_score');
                    nanIdx = any(isnan(cor),2);
                    cor(isnan(cor(:))) = 0; % set 0 if nan
                    [sorted,ord] = sort(abs(cor),'descend');
                    useVoxCand = sort(ord(1:min([p.nSelectVoxels,nVox]),:),1);
                    [uniUseVoxCand] = unique(useVoxCand(:));
                    uniUseVoxCand(ismember(uniUseVoxCand,find(nanIdx))) = [];
                    train_brain_use = train_brain(uniUseVoxCand,:)';
                    test_brain_use = test_brain(uniUseVoxCand,:);
                    
                    % add bias term
                    Xtr_1 = addBias(train_brain_use);
                    % calcurate weights
                    X1X1 = Xtr_1'*Xtr_1;
                    X1T = Xtr_1'*train_score';
                    for i = 1:nFeature
                        useVoxIdx = find(ismember(uniUseVoxCand,useVoxCand(:,i)));
                        for pitr = 1:p.nparamLogSearch
                            Lam = eye(length(useVoxIdx)+1)*p.lambda(pitr);
                            XXL = X1X1([useVoxIdx;end],[useVoxIdx;end])+Lam;
                            W = XXL\X1T([useVoxIdx;end],i);
                            pred_tmp = addBias(test_brain_use(useVoxIdx,:)')*W;
                            res_tmp.tepred_score{pitr}{ixx}(:,i) = pred_tmp.*nparm.xnorm(i) + nparm.xmean(i);
                        end
                    end
                end
                
                
                % decoding summary
                for pitr = 1:p.nparamLogSearch
                    tepred_score = merge(res_tmp.tepred_score{pitr},1); % merge prediction
                    scoreacc = fcorrdiag(L.feat(inds_all,:),tepred_score); % evaluate correlations between true and pred
                    res.scoreacc{pitr} = single(scoreacc); % reserve prediction accuracy
                    tepred_score = merge(res_tmp.tepred_score{pitr},1); % merge prediction
                    pred.prediction{pitr} = single(tepred_score); %
                end
                clear res_tmp res_ tepred*
                
                % display accuracy
                acc = zeros(p.nparamLogSearch,nFeature);
                for pitr = 1:p.nparamLogSearch
                    acc(pitr,:) = cellmean(res.scoreacc_nest(pitr,:));
                end
                [cvmx,mxind] = max(acc,[],1);
                for i = 1:nFeature
                    fprintf('score dec acc.[%s]\n',L.featname{i})
                    fprintf('r=%.3f[innerCV][lambda=%.3f]\n',cvmx(i),p.lambda(mxind(i)))
                    fprintf('r=%.3f[outerCV][lambda=%.3f]\n',res.scoreacc{mxind(i)}(i),p.lambda(mxind(i)))
                end
                
                % rearrange cv prediction accuracy
                clear cv
                cv.cvScoreAccEach_feat_bestVal = zeros(nFeature,1);
                cv.cvScoreAccEach_feat_bestIdx = zeros(nFeature,1);
                cv.cvScoreAccEach_feat = cell(nFeature,1);
                for unitr = 1:nFeature
                    for folditr = 1:p.nCVfolds
                        for paramitr = 1:p.nparamLogSearch
                            cv.cvScoreAccEach_feat{unitr}(paramitr,folditr) = res.scoreacc_nest{paramitr,folditr}(unitr);
                            cv.cvScoreAccEach_feat{unitr}(paramitr,folditr) = res.scoreacc_nest{paramitr,folditr}(unitr);
                        end
                    end
                    [cv.cvScoreAccEach_feat_bestVal(unitr),cv.cvScoreAccEach_feat_bestIdx(unitr)] = ...
                        max(mean(cv.cvScoreAccEach_feat{unitr},2));
                end
                
                % get best prediction
                cntte = 0;
                pred_best_cv_tmp = cell(p.nCVfolds,nFeature);
                for ixx = 1:p.nCVfolds
                    % get sample index
                    inds_te = find(ismember(metainf.Run,p.run2FoldAssignIdx(1,ixx):p.run2FoldAssignIdx(2,ixx))); % get test index
                    label_index = metainf.Label(inds_te);
                    inds_te = inds_te(~ismember(label_index,p.dupidx));
                    nsample_te = length(inds_te);
                    
                    % select best prediction based on each cv performance
                    for unitr = 1:nFeature
                        [val,mxind] = max(cv.cvScoreAccEach_feat{unitr}(:,folditr));
                        pred_best_cv_tmp{ixx,unitr} = double(pred.prediction{mxind}((1:nsample_te)+cntte,unitr));
                    end
                    cntte = cntte + nsample_te;
                end
                pred_best = zeros(length(inds_all),nFeature);
                for unitr = 1:nFeature
                    pred_best(:,unitr) = merge(pred_best_cv_tmp(:,unitr));
                end
                
                % identification
                clear r_best iden
                nanIdx = any(isnan(pred_best),1)+any(isnan(L.feat(inds_all,:)),1);
                ac = fcorrdiag(pred_best(:,~nanIdx),L.feat(inds_all,~nanIdx));
                fc = fcorr(pred_best(:,~nanIdx)',L.feat(inds_all,~nanIdx)');
                [iden.cr,iden.ci,iden.wins] = pwidentification(fc,1:size(pred_best,1));
                iden.cr = single(iden.cr);
                iden.ci = single(iden.ci);
                fprintf('cr(raw) = %.1f%%\n',mean(iden.cr)*100)
                
                r_best.pred_best = single(pred_best);
                r_best.all_best = ac;
                r_best.all_pattern_best = single(diag(fc));
                
                fprintf('%s\n',saveFname)
                save(saveFname,'cv','r_best','iden','-v7.3')
                tims
            end
        end
    end
end

%%
end % end function
