function draw_figure6CDE_S6(p)
%
% This code is for drawing figure 6CDE and S6
% Before drawing, clustering analysis (perform_clustering_analysis.m) must be conducted.
%
%
%% settings
testSubjectTypes = {p.sbjID{:},'randctrl'};
nclust_all = [50,27,15]; % numbers of kmeans clustering
nTopSampleProp = 5; % percentage of assigning samples
nSbj = length(p.sbjID);

cat = load([p.featdir,'category.mat']);
dim = load([p.featdir,'dimension.mat']);

% prepare scores
cat2734Score = cat.L.feat(:,[1:11,13,15:17,19:20,22:25,27:32]);
dimScore = dim.L.feat-5;

ncat2734 = size(cat2734Score,2);
ndim = size(dimScore,2);

% category and dimension names
cattxt = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','?'};
pdimtxt = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n'};
ndimtxt = {'o','p','q','r','s','t','u','v','w','x','y','z','?','+'};

%% load results
for nclitr = 1:length(nclust_all)
    nclust = nclust_all(nclitr);
    
    cnts = 0;
    for sbjitr = 1:length(testSubjectTypes)
        testSubjectType = testSubjectTypes{sbjitr};
        
        suffix_full = [testSubjectType,'/',p.suffix_enc,'/'];
        
        dataFname = sprintf('%s/kmeans_clustering/%s/kmeans_%dclust_top%d_res.mat',p.savdir,suffix_full,nclust,nTopSampleProp);
        setdir(fileparts(dataFname));
        
        switch testSubjectType
            case 'randctrl'
                tmp = load(dataFname,'H');
                
                % entorpy calculation from control
                nRand = size(tmp.H.cat,2);
                [sortedc2734] = sort(tmp.H.cat2734,2,'ascend');
                [sortedpd] = sort(tmp.H.dim_posi,2,'ascend');
                [sortednd] = sort(tmp.H.dim_nega,2,'ascend');
                pvalc2734 = 0.01/ncat2734/nSbj;
                pvald = 0.01/ndim/nSbj;
                
                D_all.ctrl.th_cat2734(nclitr) = mean(sortedc2734(:,round(nRand*pvalc2734)));
                D_all.ctrl.th_pdim(nclitr) = mean(sortedpd(:,round(nRand*pvald)));
                D_all.ctrl.th_ndim(nclitr) = mean(sortednd(:,round(nRand*pvald)));
                
            otherwise
                cnts = cnts+1;
                tmp = load(dataFname,'D','H','clustlabel');
                D_all.clustlabel{cnts,nclitr} = tmp.clustlabel;
                D_all.cat2734{cnts,nclitr} = tmp.D.cat2734;
                D_all.dim_posi{cnts,nclitr} = tmp.D.dim_posi;
                D_all.dim_nega{cnts,nclitr} = tmp.D.dim_nega;
                H_all.cat2734{cnts,nclitr} = tmp.H.cat2734;
                H_all.dim_posi{cnts,nclitr} = tmp.H.dim_posi;
                H_all.dim_nega{cnts,nclitr} = tmp.H.dim_nega;
        end
    end
end


%% draw figure 6C and S6A
% figure settings
fsize = 3;
r = 4;
c = 6;
[r,c,o] = setrc2(r*c,'ltr',[r,c],[0]);

for nclitr = find(ismember(nclust_all,27))%1:length(nclust_all)
    nclust = nclust_all(nclitr);
    h = ffigure;
    colormap bw
    cnt = 0;
    for sbjitr = 1:nSbj
        % sort while sorting both row/col
        [s,cord,rord] = sortMatDiag(D_all.cat2734{sbjitr,nclitr}',1);
        [s,cordpd,rordpd] = sortMatDiag(D_all.dim_posi{sbjitr,nclitr}',1);
        [s,cordnd,rordnd] = sortMatDiag(D_all.dim_nega{sbjitr,nclitr}',1);
        mat1 = D_all.cat2734{sbjitr,nclitr}(rord,cord);
        mat2 = D_all.dim_posi{sbjitr,nclitr}(rordpd,cordpd);
        mat3 = D_all.dim_nega{sbjitr,nclitr}(rordnd,cordnd);
        
        % category
        cnt = cnt+1;
        subplottight(r,c,o(cnt),0.2);
        imagesc(mat1',[0,30])
        set(gca,'FontSize',fsize)
        hh = colorbar;
        set(hh,'FontSize',fsize)
        axis image
        axname(cattxt(cord),2)
        ylim([0,size(mat1,2)]+0.5)
        xlabel('Cluster index')
        title(sprintf('%s:%s',p.sbjID{sbjitr},'category'))
        ffine(h)
        
        % posi dim
        cnt = cnt+1;
        subplottight(r,c,o(cnt),0.2);
        imagesc(mat2',[0,30])
        set(gca,'FontSize',fsize)
        hh = colorbar;
        set(hh,'FontSize',fsize)
        axis image
        axname(pdimtxt(cordpd),2)
        ylim([0,size(mat2,2)]+0.5)
        xlabel('Cluster index')
        title(sprintf('%s:dimPosi',p.sbjID{sbjitr}))
        ffine(h)
        
        % negadim
        cnt = cnt+1;
        subplottight(r,c,o(cnt),0.2);
        imagesc(mat3',[0,30])
        set(gca,'FontSize',fsize)
        hh = colorbar;
        set(hh,'FontSize',fsize)
        axis image
        axname(ndimtxt(cordnd),2)
        ylim([0,size(mat3,2)]+0.5)
        xlabel('Cluster index')
        title(sprintf('%s:dimNega',p.sbjID{sbjitr}))
        ffine(h)
        
    end
    suptitle(sprintf('2-D histogram of decoding accuracy for individual videos with %d kmeans cluster',nclust));
    fname = [p.figdir,'figure6C_S6A_',num2str(nclust),'cluster.pdf'];
    pathtitle(fname)
    fprintf('%s\n',savprint(h,fname));
end

%% draw figure 6D

% figure settings
fsize = 8;
r = 5;
c = 4;
[r,c,o] = setrc2(r*c,'ltd',[r,c],[1,0]);

colnd = [23,125,255];

for nclitr = 1:length(nclust_all)
    nclust = nclust_all(nclitr);
    h = ffigure;
    cnt = 0;
    
    sDcat = cell(nSbj,1);
    sDpdim = cell(nSbj,1);
    sDndim = cell(nSbj,1);
    for sbjitr = 1:nSbj
        sDcat{sbjitr} = sort(D_all.cat2734{sbjitr,nclitr},1,'descend');
        sDpdim{sbjitr} = sort(D_all.dim_posi{sbjitr,nclitr},1,'descend');
        sDndim{sbjitr} = sort(D_all.dim_nega{sbjitr,nclitr},1,'descend');
    end
    
    % individual results
    for sbjitr = 1:nSbj
        % sorted histogram
        cnt = cnt+1;
        subplottight(r,c,o(cnt),0.2);
        hh = plot(sDcat{sbjitr},'Color',[255,69,102]/255);
        set(gca,'FontSize',fsize)
        AnnotationOff(hh,2:length(hh));
        hold on
        hh = plot(sDpdim{sbjitr},'Color',[23,125,255]/255);
        AnnotationOff(hh,2:length(hh));
        hold on
        plot(sDndim{sbjitr},'Color',colnd/255)
        
        axname([1,5:5:25],1,[1,5:5:25])
        xlim([0,nclust+1])
        ylim([0,45]);
        hh = hline(0:5:45,'-w');
        
        set(hh,'Color',[1,1,1]*0.7)
        ylabel('Frequency')
        xlabel('Sorted cluster index')
        legend({'category','dimPosi','dimNega'})
        title(sprintf('%s',p.sbjID{sbjitr}))
        ffine(h)
    end
    
    % averaged results
    cnt = cnt+1;
    subplottight(r,c,o(cnt),0.2);
    hh = plot(cellmean(sDcat),'Color',[255,69,102]/255);
    set(gca,'FontSize',fsize)
    AnnotationOff(hh,2:length(hh));
    hold on
    plot(cellmean(sDpdim),'Color',[23,125,255]/255)
    hold on
    
    plot(cellmean(sDndim),'Color',colnd/255)
    xlim([0,nclust+1])
    ylim([0,45]);
    hh = hline(0:5:45,'-w');
    
    set(hh,'Color',[1,1,1]*0.7)
    ylabel('Frequency')
    xlabel('Sorted cluster index')
    legend({'category','dimPosi','dimNega'})
    title(sprintf('%s:','Average'))
    ffine(h)
    
    
    suptitle(sprintf('1-D histogram of decoding accuracy for individual videos with %d kmeans cluster',nclust));
    fname = [p.figdir,'figure6D_S6BD_',num2str(nclust),'cluster.pdf'];
    pathtitle(fname)
    fprintf('%s\n',savprint(h,fname));
end

%% draw figure 6E

% figure settings
fsize = 8;
fsizeletter = 5;
r = 5;
c = 4;
[r,c,o] = setrc2(r*c,'ltd',[r,c],[1,0]);
msize = 1;


for nclitr = 1:length(nclust_all)
    nclust = nclust_all(nclitr);
    h = ffigure;
    cnt = 0;
    
    thcat = D_all.ctrl.th_cat2734(nclitr);
    th_pdim = D_all.ctrl.th_pdim(nclitr);
    th_ndim = D_all.ctrl.th_ndim(nclitr);
    
    for sbjitr = 1:nSbj
        cnt = cnt+1;
        subplottight(r,c,o(cnt),0.2);
        sinaplot({thcat,th_pdim,th_ndim},'MarkerType','none','mc','r','ci',0);
        hold on
        hh = sinaplot({H_all.cat2734{sbjitr,nclitr},H_all.dim_posi{sbjitr,nclitr},H_all.dim_nega{sbjitr,nclitr}},...
            'MarkerType','none','edgecolor',[1,1,1],'MarkerSize',msize,'mc',[],'ci',0);
        set(gca,'FontSize',fsize)
        
        % color only for significan ROIs
        xd = get(hh,'XData');
        yd = get(hh,'YData');
        
        % category
        col = cmap4('ck27');
        for i = 1:length(xd{1})
            text(xd{1}(i),yd{1}(i),cattxt{i},'Color',col{i},'FontSize',fsizeletter,'HorizontalAlignment','center');
        end
        % posi dim
        col = cmap4('ck28c');
        col = col(1:14,:);
        for i = 1:length(xd{2})
            text(xd{2}(i),yd{2}(i),pdimtxt{i},'Color',col{i},'FontSize',fsizeletter,'HorizontalAlignment','center');
        end
        % nega dim
        col = cmap4('ck28c');
        col = col(15:28,:);
        for i = 1:length(xd{3})
            text(xd{3}(i),yd{3}(i),ndimtxt{i},'Color',col{i},'FontSize',fsizeletter,'HorizontalAlignment','center');
        end
        switch nclust
            case 50
                ylim([3.5,5.5])
                hh = hline(3.5:0.5:5.5,'-w');
                axname(3.5:0.5:5.5,2,3.5:0.5:5.5)
            case {15}
                ylim([2,4])
                hh = hline(2:0.4:4,'-w');
                axname(2:0.4:4,2,2:0.4:4)
            case {27}
                ylim([2.5,5])
                hh = hline(2.5:0.5:5,'-w');
                axname(2.5:0.5:5,2,2.5:0.5:5)
        end
        ylabel('Entropy')
        set(hh,'Color',[1,1,1]*0.3)
        axname({'category','dimPosi','dimPNega'})
        xticklabel_rotate([],45)
        ppd = ranksum(H_all.cat2734{sbjitr,nclitr},H_all.dim_posi{sbjitr,nclitr});
        pnd = ranksum(H_all.cat2734{sbjitr,nclitr},H_all.dim_nega{sbjitr,nclitr});
        pd = ranksum(H_all.cat2734{sbjitr,nclitr},[H_all.dim_posi{sbjitr,nclitr};H_all.dim_nega{sbjitr,nclitr}]);
        title(sprintf('%s:%s:[ppd=%.4f,pnd=%.4f,pd=%.4f]',p.sbjID{sbjitr},'category',ppd,pnd,pd))
        ffine(h)
    end
    
    % Average
    cnt = cnt+1;
    subplottight(r,c,o(cnt),0.2);
    sinaplot({thcat,th_pdim,th_ndim},'MarkerType','none','mc','r','ci',0);
    hold on
    hh = sinaplot({cellmean(H_all.cat2734(:,nclitr)),cellmean(H_all.dim_posi(:,nclitr)),cellmean(H_all.dim_nega(:,nclitr))},...
        'MarkerType','none','edgecolor',[1,1,1],'MarkerSize',msize,'mc',[],'ci',0);
    set(gca,'FontSize',fsize)
    
    % color only for significan ROIs
    xd = get(hh,'XData');
    yd = get(hh,'YData');
    
    % category
    col = cmap4('ck27');
    for i = 1:length(xd{1})
        text(xd{1}(i),yd{1}(i),cattxt{i},'Color',col{i},'FontSize',fsizeletter,'HorizontalAlignment','center');
    end
    % posi dim
    col = cmap4('ck28c');
    col = col(1:14,:);
    for i = 1:length(xd{2})
        text(xd{2}(i),yd{2}(i),pdimtxt{i},'Color',col{i},'FontSize',fsizeletter,'HorizontalAlignment','center');
    end
    % nega dim
    col = cmap4('ck28c');
    col = col(15:28,:);
    for i = 1:length(xd{3})
        text(xd{3}(i),yd{3}(i),ndimtxt{i},'Color',col{i},'FontSize',fsizeletter,'HorizontalAlignment','center');
    end
    
    switch nclust
        case 50
            ylim([3.5,5.5])
            hh = hline(3.5:0.5:5.5,'-w');
            axname(3.5:0.5:5.5,2,3.5:0.5:5.5)
        case {15}
            ylim([2,4])
            hh = hline(2:0.4:4,'-w');
            axname(2:0.4:4,2,2:0.4:4)
        case {27}
            ylim([2.5,5])
            hh = hline(2.5:0.5:5,'-w');
            axname(2.5:0.5:5,2,2.5:0.5:5)
    end
    
    ylabel('Entropy')
    set(hh,'Color',[1,1,1]*0.3)
    axname({'category','dimPosi','dimPNega'})
    xticklabel_rotate([],45)
    ppd = ranksum(cellmean(H_all.cat2734),cellmean(H_all.dim_posi));
    pnd = ranksum(cellmean(H_all.cat2734),cellmean(H_all.dim_nega));
    pd = ranksum(cellmean(H_all.cat2734),[cellmean(H_all.dim_posi);cellmean(H_all.dim_nega)]);
    title(sprintf('%s:%s:[ppd=%.4f,pnd=%.4f,pd=%.4f]','Average','category',ppd,pnd,pd))
    ffine(h)
    
    
    suptitle(sprintf('Entropy distributions with %d kmeans cluster',nclust));
    fname = [p.figdir,'figure6E_S6CE_',num2str(nclust),'cluster.pdf'];
    pathtitle(fname)
    fprintf('%s\n',savprint(h,fname));
end


%%
close all
