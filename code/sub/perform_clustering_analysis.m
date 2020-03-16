function perform_clustering_analysis(p,encRes)
%
% - This code is written for performing clustering analysis
%   based on encoding accuracy
%
% called from emotion2020_analysis_BATCH
%

%% settings
testSubjectTypes = {p.sbjID{:},'randctrl'};

nclust_all = [50,27,15]; % numbers of kmeans clustering
nTopSampleProp = 5; % percentage of assigning samples
cvth = 0.124; % threshold for voxel selection
selectVoxProp = 100; % proportion of voxels used for voxel selection

nSample = 2181;
nRand = 100000; % # of shuffling computations


cat = load([p.featdir,'category.mat']);
dim = load([p.featdir,'dimension.mat']);


%% perform kmeans clustering
for sbjitr = 1:length(testSubjectTypes)
    testSubjectType = testSubjectTypes{sbjitr};
    
    suffix_full = [testSubjectType,'/',p.suffix_enc,'/'];
    
    for nclitr = randsample(1:length(nclust_all),length(nclust_all))%1:length(nclust_all)
        nclust = nclust_all(nclitr);
        
        
        saveFname = sprintf('%s/kmeans_clustering/%s/kmeans_%dclust_top%d_res.mat',p.savdir,suffix_full,nclust,nTopSampleProp);
        saveFnameChk = sprintf('%s/kmeans_clustering/%s/kmeans_%dclust_top%d_log.txt',p.logdir,suffix_full,nclust,nTopSampleProp);
        setdir(fileparts(saveFname));
        setdir(fileparts(saveFnameChk));
        
        if exist(saveFnameChk,'file')
            fprintf('%s\n',saveFname)
        else
            fprintf('Start:%s\n',saveFnameChk)
            tmp = [];
            save(saveFnameChk,'tmp','-ascii')
            
            % subject
            switch testSubjectType
                case 'randctrl'
                    
                    % prepare scores
                    catScore = cat.L.feat;
                    dimScore = dim.L.feat-5;
                    cat2734Score = cat.L.feat(:,[1:11,13,15:17,19:20,22:25,27:32]);
                    
                    % remove duplicates
                    catScore(p.dupidx,:) = [];
                    dimScore(p.dupidx,:) = [];
                    cat2734Score(p.dupidx,:) = [];
                    
                    [sorted,ordc] = sort(catScore,'descend');
                    [sorted,ordc2734] = sort(cat2734Score,'descend');
                    [sorted,ordd] = sort(dimScore,'descend');
                    
                    n = round(nSample/100*nTopSampleProp);
                    
                    clear H D
                    
                    ncat = size(catScore,2);
                    ncat2734 = size(cat2734Score,2);
                    ndim = size(dimScore,2);
                    
                    D.cat = cell(nRand,1);
                    D.cat2734 = cell(nRand,1);
                    D.dim_posi = cell(nRand,1);
                    D.dim_nega = cell(nRand,1);
                    H.cat = zeros(ncat,nRand);
                    H.cat2734 = zeros(ncat2734,nRand);
                    H.dim_posi = zeros(ndim,nRand);
                    H.dim_nega = zeros(ndim,nRand);
                    
                    for randitr = 1:nRand
                        rand('seed',randitr)
                        clustlabel = randsample(1:nclust,nSample,1)';
                        
                        
                        D.cat{randitr} = zeros(nclust,ncat);
                        D.cat2734{randitr} = zeros(nclust,ncat2734);
                        D.dim_posi{randitr} = zeros(nclust,ndim);
                        D.dim_nega{randitr} = zeros(nclust,ndim);
                        
                        for ix = 1:ncat
                            D.cat{randitr}(:,ix) = histc(clustlabel(ordc(1:n,ix)),1:nclust)';
                        end
                        for ix = 1:ncat2734
                            D.cat2734{randitr}(:,ix) = histc(clustlabel(ordc2734(1:n,ix)),1:nclust)';
                        end
                        for ix = 1:ndim
                            D.dim_posi{randitr}(:,ix) = histc(clustlabel(ordd(1:n,ix)),1:nclust)';
                            D.dim_nega{randitr}(:,ix) = histc(clustlabel(ordd(end-n+1:end,ix)),1:nclust)';
                        end
                        
                        for ix = 1:ncat
                            distrib = histc(clustlabel(ordc(1:n,ix)),1:nclust)';
                            H.cat(ix,randitr) = entropy0(distrib./nclust);
                        end
                        for ix = 1:ncat2734
                            distrib = histc(clustlabel(ordc2734(1:n,ix)),1:nclust)';
                            H.cat2734(ix,randitr) = entropy0(distrib./nclust);
                        end
                        for ix = 1:ndim
                            pdistrib = histc(clustlabel(ordd(1:n,ix)),1:nclust)';
                            H.dim_posi(ix,randitr) = entropy0(pdistrib./nclust);
                            
                            ndistrib = histc(clustlabel(ordd(end-n+1:end,ix)),1:nclust)';
                            H.dim_nega(ix,randitr) = entropy0(ndistrib./nclust);
                        end
                        
                        if mod(randitr,nRand/100) == 0
                            fprintf('%d/%d\n',randitr,nRand)
                            s = sort(asvector(H.cat(:,1:randitr)));
                            fprintf('cat :[p=0.01],H = %.2f,[min],H = %.5f\n',s(ceil(length(s)/100)),min(s))
                            s = sort(asvector(H.cat2734(:,1:randitr)));
                            fprintf('cat2734:[p=0.01],H = %.2f,[min],H = %.5f\n',s(ceil(length(s)/100)),min(s))
                            s = sort(asvector(H.dim_posi(:,1:randitr)));
                            fprintf('pdim:[p=0.01],H = %.2f,[min],H = %.5f\n',s(ceil(length(s)/100)),min(s))
                            s = sort(asvector(H.dim_nega(:,1:randitr)));
                            fprintf('ndim:[p=0.01],H = %.2f,[min],H = %.5f\n',s(ceil(length(s)/100)),min(s))
                        end
                    end
                    fprintf('%s\n',saveFname)
                    save(saveFname,'D','H','-v7.3')
                    
                otherwise
                    
                    fprintf('Load data ...\n')
                    dpath = sprintf('%s%s/rois/%s_WholeBrain.mat',p.fmridir,p.sbjID{sbjitr},p.sbjID{sbjitr});
                    load(dpath,'braindat','metainf');
                    braindat = zscore(braindat);
                    % brain data info.
                    nSample = size(braindat,1);
                    n = round(nSample/100*nTopSampleProp);
                    
                    % voxel selection
                    mDimCV = mean(encRes.dimension.pred_acc_cv{sbjitr},1);
                    mCatCV = mean(encRes.category.pred_acc_cv{sbjitr},1);
                    mSemCV = mean(encRes.semantic.pred_acc_cv{sbjitr},1);
                    mVisCV = mean(encRes.vision.pred_acc_cv{sbjitr},1);
                    
                    % select voxels predicted by emotion models with significant accuracy
                    [mx,mxind] = max([mCatCV;mDimCV;mSemCV;mVisCV],[],1);
%                     voxSelIndx = (mCatCV > cvth|mDimCV > cvth) & (mxind == 1 | mxind == 2);
                    voxSelIndx = (mCatCV > cvth|mDimCV > cvth);% & (mxind == 1 | mxind == 2);
%                     [sorted,ord] = sort(mCatCV,'descend');
%                     topVoxIdx = ord(1:sum(voxSelIndx)*selectVoxProp/100);
%                     voxSelIndx = false(size(voxSelIndx));
%                     voxSelIndx(topVoxIdx) = true;
                    
                    % sort labels
                    label_index = metainf.Label;
                    [sorted,ord] = sort(label_index);
                    b_vox = braindat(ord,voxSelIndx);
                    b_vox(p.dupidx,:) = [];
                    
                    corrdist_vox = 1-fcorr(b_vox',b_vox');
                    
                    % prepare scores
                    catScore = cat.L.feat;
                    dimScore = dim.L.feat-5;
                    cat2734Score = cat.L.feat(:,[1:11,13,15:17,19:20,22:25,27:32]);
                    
                    % remove duplicates
                    catScore(p.dupidx,:) = [];
                    dimScore(p.dupidx,:) = [];
                    cat2734Score(p.dupidx,:) = [];
                    
                    [sorted,ordc] = sort(catScore,'descend');
                    [sorted,ordc2734] = sort(cat2734Score,'descend');
                    [sorted,ordd] = sort(dimScore,'descend');
                    
                    ncat = size(catScore,2);
                    ncat2734 = size(cat2734Score,2);
                    ndim = size(dimScore,2);
                    
                    opts = statset('Display','final');
                    [clustlabel,C] = kmeans(b_vox,nclust,'Distance','correlation',...
                        'Replicates',10,'Options',opts);
                    
                    
                    clear H D
                    D.cat = zeros(nclust,ncat);
                    D.cat2734 = zeros(nclust,ncat2734);
                    D.dim_posi = zeros(nclust,ndim);
                    D.dim_nega = zeros(nclust,ndim);
                    
                    H.cat = zeros(ncat,1);
                    H.cat2734 = zeros(ncat2734,1);
                    H.dim_posi = zeros(ndim,1);
                    H.dim_nega = zeros(ndim,1);
                    
                    
                    for ix = 1:ncat
                        D.cat(:,ix) = histc(clustlabel(ordc(1:n,ix)),1:nclust)';
                    end
                    for ix = 1:ncat2734
                        D.cat2734(:,ix) = histc(clustlabel(ordc2734(1:n,ix)),1:nclust)';
                    end
                    for ix = 1:ndim
                        D.dim_posi(:,ix) = histc(clustlabel(ordd(1:n,ix)),1:nclust)';
                        D.dim_nega(:,ix) = histc(clustlabel(ordd(end-n+1:end,ix)),1:nclust)';
                    end
                    
                    for ix = 1:ncat
                        distrib = histc(clustlabel(ordc(1:n,ix)),1:nclust)';
                        H.cat(ix) = entropy0(distrib./nclust);
                    end
                    for ix = 1:ncat2734
                        distrib = histc(clustlabel(ordc2734(1:n,ix)),1:nclust)';
                        H.cat2734(ix) = entropy0(distrib./nclust);
                    end
                    for ix = 1:ndim
                        pdistrib = histc(clustlabel(ordd(1:n,ix)),1:nclust)';
                        H.dim_posi(ix) = entropy0(pdistrib./nclust);
                        
                        ndistrib = histc(clustlabel(ordd(end-n+1:end,ix)),1:nclust)';
                        H.dim_nega(ix) = entropy0(ndistrib./nclust);
                    end
                    
                    fprintf('%s\n',saveFname)
                    save(saveFname,'D','H','b_vox','corrdist_vox','label_index','clustlabel','-v7.3')
                    
            end
        end
    end
end

%%

