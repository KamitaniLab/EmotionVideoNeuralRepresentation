function perform_encoding_analysis(p)
%
% - This code is written for performing whole brain encoding analysis from specified features.
% - Whole brain data was parsed into "nparse" groups for parallel computations.
% - A whole brain analysis for one subject and for one feature type takes about 2-3 hours with single cpu.
%
%
% called from emotion2020_analysis_BATCH
%
%% initialization
warning off

%% start analysis
for scoritr = 1:length(p.scoreTypes)
    scoreType = p.scoreTypes{scoritr};
    
    for sbjitr = 1:length(p.sbjID)
        sbjID = p.sbjID{sbjitr};
        
        % set suffix
        suffix_full = [sbjID,'/',p.suffix_enc];
        
        %% Save info. for summary results
        saveFnameChk         = sprintf('%s/%s/%s/_log.txt',p.logdir,scoreType,suffix_full); % log files
        saveFname            = sprintf('%s/%s/%s/res.mat',p.savdir,scoreType,suffix_full); % res files
        setdir(fileparts(saveFnameChk));
        setdir(fileparts(saveFname));
        
        %% train-test analyses
        tmp = []; % empty variable for log file
        fCheck = zeros(p.nparse,1); % check variable
        cnt = 0; % count variable
        
        if p.checkModeRes
            chkfile = saveFname;
        else
            chkfile = saveFnameChk;
        end
        if exist(chkfile,'file')
            fCheck = zeros(p.nparse,1); % check whther the subset of analysis was finished
            if p.checkChkfile
                if exist(saveFname,'file')
                    fCheck = ones(p.nparse,1);
                else
                    for itr = 1:p.nparse
                        saveFnameChkx = sprintf('%s/%s/%s/_p%02d_log.txt',p.logdir,scoreType,suffix_full,itr);
                        saveFnamex    = sprintf('%s/%s/%s/_p%02d.mat',p.savdir,scoreType,suffix_full,itr);
                        if exist(saveFnamex,'file')
                            fCheck(itr) = 1;
                        else
                            if p.del
                                fprintf('Delete%s\n',saveFnameChkx)
                                delete(saveFnameChkx)
                                fprintf('Delete%s\n',saveFnameChk)
                                delete(saveFnameChk)
                            end
                            fCheck(itr) = 0;
                        end
                    end
                end
            end
        else
            tic
            
            % load brain data and get params
            fprintf('Load data(%s)...\n',sbjID)
            dpath = sprintf('%s%s/rois/%s_WholeBrain.mat',p.fmridir,sbjID,sbjID);
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
            
            
            
            %% loop for small subsets
            for itr = randsample(1:p.nparse,p.nparse) % 1:p.nparse is OK
                tic
                cnt = cnt+1;
                saveFnameChkx = sprintf('%s/%s/%s/_p%02d_log.txt',p.logdir,scoreType,suffix_full,itr);
                saveFnamex    = sprintf('%s/%s/%s/_p%02d.mat',p.savdir,scoreType,suffix_full,itr);
                if p.checkModeRes
                    chkfilex = saveFnamex;
                else
                    chkfilex = saveFnameChkx;
                end
                if exist(chkfilex,'file')
                    if exist(saveFnamex,'file')
                        fCheck(itr) = 1;
                    else
                        fCheck(itr) = 0;
                        if p.del
                            fprintf('Delete%s\n',saveFnameChkx)
                            delete(saveFnameChkx)
                        end
                    end
                elseif ~p.del
                    fprintf('Start:%s\n',saveFnameChkx)
                    fprintf('CV %s %02d parse[%d/%d][%s]:\n',suffix_full,itr,cnt,p.nparse,scoreType)
                    save(saveFnameChkx,'tmp', '-ascii') % save log file
                    
                    % voxel group selection
                    if p.nparse>1
                        n = ceil(nVox/p.nparse);
                        if nVox/p.nparse < 1
                            idx = (1:nVox)+(itr-1)*nVox;
                        else
                            idx = (1:n)+(itr-1)*n;
                        end
                        idx = idx(idx<=nVox);
                    else
                        idx = 1:nVox;
                    end
                    
                    % set containers for a voxel group
                    clear res pred pred_nest
                    res.bdacc.all = cell(p.nparamLogSearch,1);
                    pred.prediction = cell(p.nparamLogSearch,1);
                    res_tmp.tepred_brain.all = cell(p.nparamLogSearch,1);
                    
                    if length(idx) < 1
                        fprintf('# of assigned voxels is 0.\n')
                    else
                        
                        
                        %% analysis loop
                        for ixx = 1:p.nCVfolds
                            % get sample index
                            inds_te = find(ismember(metainf.Run,p.run2FoldAssignIdx(1,ixx):p.run2FoldAssignIdx(2,ixx))); % get test index
                            label_index = metainf.Label(inds_te);
                            inds_te = inds_te(~ismember(label_index,p.dupidx));
                            inds_tr = inds_all(~ismember(inds_all,inds_te));
                            
                            clear D_tr D_te parm res_tmp_nest
                            % separate training and test data
                            D_tr.brain = braindat(inds_tr,idx);
                            D_te.brain = braindat(inds_te,idx);
                            D_tr.score = L.feat(inds_tr,:);
                            D_te.score = L.feat(inds_te,:);
                            
                            % nested CV ===============================
                            % nested CV was performed within the Training data
                            trainRunIdx = ~ismember(1:p.nCVfolds,ixx);
                            nestcvIdx = [p.run2FoldAssignIdx(1,trainRunIdx);p.run2FoldAssignIdx(2,trainRunIdx)];
                            nCVfolds_nest = size(nestcvIdx,2);
                            Run_train = metainf.Run(inds_tr);
                            inds_all_nest = 1:length(inds_tr);
                            
                            % container for each voxel.
                            res_tmp_nest.tepred_brain.all = cell(p.nparamLogSearch,1);
                            
                            for ixx2 = 1:nCVfolds_nest
                                clear D_tr_nest D_te_nest
                                % separate training and test data
                                inds_te_nest = find(ismember(Run_train,nestcvIdx(1,ixx2):nestcvIdx(2,ixx2))); % get test index
                                inds_tr_nest = find(~ismember(inds_all_nest,inds_te_nest)); % get traiing index
                                D_tr_nest.brain = D_tr.brain(inds_tr_nest,:);
                                D_te_nest.brain = D_tr.brain(inds_te_nest,:);
                                D_tr_nest.score = D_tr.score(inds_tr_nest,:);
                                D_te_nest.score = D_tr.score(inds_te_nest,:);
                                
                                % normalization
                                % score normalization
                                parm.data_norm = 1;
                                [train_score,nparm] = normalize_data(D_tr_nest.score',parm.data_norm);
                                parm.xmean = nparm.xmean;
                                parm.xnorm = nparm.xnorm;
                                test_score = normalize_data(D_te_nest.score',parm.data_norm,parm);
                                
                                % brain normalization
                                [train_brain,nparm] = normalize_data(double(D_tr_nest.brain)',parm.data_norm);
                                parm_all_nest.ymean=nparm.xmean;
                                parm_all_nest.ynorm=nparm.xnorm;
                                
                                % perform l2 regression
                                % training & test
                                nsamp = size(test_score,2);
                                % add the intercept
                                Xtr_1 = addBias(train_score');
                                X1X1 = Xtr_1'*Xtr_1;
                                % training & test
                                for pitr = 1:p.nparamLogSearch
                                    [pred_tmp,weight] = flinreg_l2(p.lambda(pitr),X1X1,Xtr_1,train_brain',test_score');
                                    res_tmp_nest.tepred_brain.all{pitr}{ixx2} = pred_tmp'.*repmat(parm_all_nest.ynorm,1,nsamp) + repmat(parm_all_nest.ymean,1,nsamp);
                                end
                            end
                            
                            % encoding summary
                            for pitr = 1:p.nparamLogSearch
                                pred_temp = merge(res_tmp_nest.tepred_brain.all{pitr},2)';
                                bdacc = fcorrdiag(D_tr.brain,pred_temp); % evaluate correlations between true and pred
                                res.bdacc.all_nest{pitr}(ixx,:) = single(bdacc); % reserve prediction accuracy
                            end
                            clear D_tr_nest D_te_nest parm parm_all_nest train_score test_score train_brain test_brain res_tmp_nest weight Model
                            % end nested CV analyses
                            % ========================
                            
                            
                            % res-start CV analyses (not nestedCV)========================
                            % normalization
                            % score normalization
                            parm.data_norm = 1;
                            [train_score,nparm] = normalize_data(D_tr.score',parm.data_norm);
                            parm.xmean = nparm.xmean;
                            parm.xnorm = nparm.xnorm;
                            test_score = normalize_data(D_te.score',parm.data_norm,parm);
                            
                            % brain normalization
                            [train_brain,nparm] = normalize_data(double(D_tr.brain)',parm.data_norm);
                            parm_all.ymean = nparm.xmean;
                            parm_all.ynorm = nparm.xnorm;
                            
                            % perform l2 regression
                            % training & test
                            nsamp = size(test_score,2);
                            % add the intercept
                            Xtr_1 = addBias(train_score');
                            X1X1 = Xtr_1'*Xtr_1;
                            % training & test
                            for pitr = 1:p.nparamLogSearch
                                [pred_tmp,res_tmp.weights{pitr}{ixx}] = flinreg_l2(p.lambda(pitr),X1X1,Xtr_1,train_brain',test_score');
                                res_tmp.tepred_brain.all{pitr}{ixx} = pred_tmp'.*repmat(parm_all.ynorm,1,nsamp) + repmat(parm_all.ymean,1,nsamp);
                            end
                        end
                        
                        % encoding summary
                        for pitr = 1:p.nparamLogSearch
                            tepred_brain = merge(res_tmp.tepred_brain.all{pitr},2)'; % merge prediction
                            bdacc = fcorrdiag(braindat(inds_all,idx),tepred_brain); % evaluate correlations between true and pred
                            res.bdacc.all{pitr} = single(bdacc)'; % reserve prediction accuracy
                            pred.prediction{pitr} = single(tepred_brain); % reserve prediction
                        end
                    end
                    % save resutls for one voxel group
                    save(saveFnamex,'res','pred','-v7.3')
                    
                    % display cs and accuracy
                    if isempty(idx)
                        fprintf('# of assigned voxels is 0.\n')
                    else
                        try
                            fprintf('[vox:%d-%d/%d(%d/%d)]\n',idx(1),idx(end),nVox,cnt,p.nparse)
                            cvacc = mean(merge(cellfun(@mean,res.bdacc.all_nest,'UniformOutput',false)),2);
                            [mx,mxind] = max(cvacc);
                            fprintf('enc acc(innerCV): r=%.3f[lambda = %.3f]\n',nanmean(nanmean(res.bdacc.all_nest{mxind}(:))),p.lambda(mxind))
                            fprintf('enc acc(outerCV): r=%.3f[lambda = %.3f]\n',nanmean(cellmean(res.bdacc.all(mxind))),p.lambda(mxind))
                        end
                    end
                    tims
                    clear res_tmp res_ tepred*  res stats pred
                end
            end
            
            % check whether all voxles groups were finished
            for itr = 1:p.nparse
                saveFnamex    = sprintf('%s/%s/%s/_p%02d.mat',p.savdir,scoreType,suffix_full,itr);
                if exist(saveFnamex,'file')
                    fCheck(itr) = 1;
                else
                    fCheck(itr) = 0;
                end
            end
        end
        
        % aggregate all resutls from all voxel groups
        if all(fCheck)&&~exist(chkfile,'file') && ~p.del && p.integrateRes
            fprintf('Integrate %s\n',saveFname)
            save(saveFnameChk,'tmp', '-ascii')
            
            % load data and get params for all subjects
            fprintf('Load metainf (%s)...\n',sbjID)
            dpath = sprintf('%s%s/rois/%s_WholeBrain.mat',p.fmridir,sbjID,sbjID);
            load(dpath,'metainf');
            
            clear result stats preds
            result_tmp.bdacc.all = cell(p.nparamLogSearch,p.nparse);
            result_tmp.bdacc.all_nest = cell(p.nparamLogSearch,p.nparse);
            result.bdacc.all = cell(p.nparamLogSearch,1);
            result.bdacc.all_nest = cell(p.nparamLogSearch,1);
            preds = cell(p.nparamLogSearch,1);
            preds_tmp = cell(p.nparamLogSearch,p.nparse);
            
            for itr = 1:p.nparse
                if mod(itr,p.nparse/10)==0
                    fprintf('%d/%d\n',itr,p.nparse)
                    tims
                end
                saveFnameChkx = sprintf('%s/%s/%s/_p%02d_log.txt',p.logdir,scoreType,suffix_full,itr);
                saveFnamex    = sprintf('%s/%s/%s/_p%02d.mat',p.savdir,scoreType,suffix_full,itr);
                clear res stats pred
                try
                    load(saveFnamex,'res','pred');
                    
                    % keep all results
                    for pitr = 1:p.nparamLogSearch
                        if isempty(res.bdacc.all)
                            fprintf('Empty results:%d(p.lambda=%.3f)\n',itr,p.lambda(pitr))
                        else
                            preds_tmp{pitr,itr} = pred.prediction{pitr};
                            result_tmp.bdacc.all{pitr,itr} = res.bdacc.all{pitr};
                            result_tmp.bdacc.all_nest{pitr,itr} = res.bdacc.all_nest{pitr};
                        end
                    end
                catch me
                    fprintf('Something wrong for loading results\n')
                    fprintf('Delete: %s\n',saveFnameChkx)
                    fprintf('Delete: %s\n',saveFnamex)
                    %delete(saveFnameChkx)
                    %delete(saveFnamex)
                    error('Try again.')
                end
            end
            % encoding summary
            clear tmp
            % remove empty group
            for pitr = 1:p.nparamLogSearch
                tmp.bdacc.all = result_tmp.bdacc.all(pitr,~cellfun(@isempty,result_tmp.bdacc.all(pitr,:))); % remove empty
                result.bdacc.all{pitr} = [tmp.bdacc.all{:}];% merge
                
                tmp.bdacc_nest = result_tmp.bdacc.all_nest(pitr,~cellfun(@isempty,result_tmp.bdacc.all_nest(pitr,:))); % remove empty
                result.bdacc.all_nest{pitr} = [tmp.bdacc_nest{:}];% merge
                
                tmp.preds = preds_tmp(pitr,~cellfun(@isempty,preds_tmp(pitr,:))); % remove empty
                preds{pitr} = [tmp.preds{:}];% merge
            end
            clear tmp

            % cv accuracy summary
            mean_acc_nest = nanmean(merge(cellfun(@nanmean,result.bdacc.all_nest,'UniformOutput',false)),2);
            [mx,mxind] = max(mean_acc_nest);
            pred_acc_cv = result.bdacc.all_nest{mxind};

            % select best params for each CV
            mean_acc_nest = zeros(p.nparamLogSearch,p.nCVfolds);
            mxindcv = zeros(p.nCVfolds,1);
            pred_best_cv = cell(p.nCVfolds,1);
            
            cntte = 0;
            % get prediction with best reguralization parameters for each CV.
            for ixx = 1:p.nCVfolds
                inds_te = find(ismember(metainf.Run,p.run2FoldAssignIdx(1,ixx):p.run2FoldAssignIdx(2,ixx))); % get test index
                label_index = metainf.Label(inds_te);
                inds_te = inds_te(~ismember(label_index,p.dupidx));
                nsample_te = length(inds_te);
                for pitr = 1:p.nparamLogSearch
                    mean_acc_nest(pitr,ixx) = nanmean(result.bdacc.all_nest{pitr}(ixx,:));
                end
                [mx,mxindcv(ixx)] = max(mean_acc_nest(:,ixx));
                pred_best_cv{ixx} = preds{mxindcv(ixx)}((1:nsample_te)+cntte,:);
                cntte = cntte + nsample_te;
            end
            pred_best = merge(pred_best_cv,1);
            clear pred_best_cv
            
            fprintf('Identification analysis...\n')
            dpath = sprintf('%s%s/rois/%s_WholeBrain.mat',p.fmridir,sbjID,sbjID);
            load(dpath,'braindat');
            braindat = zscore(braindat);
            
            % remove duplicates
            inds_all = 1:size(braindat,1);
            label_index = metainf.Label;
            inds_all = inds_all(~ismember(label_index,p.dupidx));
            
            clear iden r_best
            nanIdx = any(isnan(pred_best),1)+any(isnan(braindat(inds_all,:)),1);
            ac = fcorrdiag(pred_best(:,~nanIdx),double(braindat(inds_all,~nanIdx)));
            r_best.pred_acc = single(ac');
            r_best.pred_acc_cv = pred_acc_cv;
            
            fc = fcorr(pred_best(:,~nanIdx)',double(braindat(inds_all,~nanIdx))');
            [r_best.iden_acc] = pwidentification(fc,1:size(pred_best,1));
            fprintf('cr = %.1f%%\n',mean(r_best.iden_acc)*100)
            
            fprintf('%s\n',saveFname)
            save(saveFname,'result','pred_best','r_best','-v7.3')
        end
    end
end

%% end function
warning on
end
