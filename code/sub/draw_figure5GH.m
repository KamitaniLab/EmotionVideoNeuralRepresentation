function draw_figure5GH(p,encRes)
%
% This code is for drawing figure 5GH
%
%
%% settings

% prepare colors
clear col
col{1}=[10, 168, 247]/255; %
col{2}=[18, 94, 240]/255; %
col{3}=[235, 1, 22]/255; %

nSbj = length(p.sbjID);

for sbjitr = 1:nSbj
    pgPath = sprintf('%s%s/pringrad/%s_pringrad_values.mat',p.fmridir,p.sbjID{sbjitr},p.sbjID{sbjitr});
    load(pgPath,'pg')
    pgMaps.(p.sbjID{sbjitr}) = pg;
end


%% draw figure 5G
compSets = {'vision','semantic','category'};

% display settings
r = 6;
c = 8;
[r,c,o] = setrc2(r*c,'ltd',[r,c],[1]);
fsize = 4;
msize = 2;%1.5;
thval = 0.111;
nPlotVoxel = 5000; % 5000 voxels were plotted in the paper

close all
h = ffigure;
cnt = 0;

accs = cell(length(compSets),nSbj);
sigs = cell(length(compSets),nSbj);
for sbjitr = 1:nSbj
    useIdx = pgMaps.(p.sbjID{sbjitr}).cortexIdx;
    for compitr = 1:length(compSets)
        accs{compitr,sbjitr} = encRes.(compSets{compitr}).pred_acc{sbjitr}(useIdx)';
        sigs{compitr,sbjitr} = accs{compitr,sbjitr} >= thval;
    end
end

for sbjitr = 1:nSbj
    acc = merge(accs(:,sbjitr),2);
    sig = any(merge(sigs(:,sbjitr),2),2);
    [mx,mxidx] = max(acc,[],2);
    
    % select only signif voxel
    acc = acc(sig,:);
    mxidx = mxidx(sig);
    nSigVox = length(acc);
    pg1 = pgMaps.(p.sbjID{sbjitr}).map{1}(sig);
    pg2 = -pgMaps.(p.sbjID{sbjitr}).map{2}(sig);
    
    % randm selection for visualization
    rand('seed',1);
    randIdx = randsample(nSigVox,nPlotVoxel);
    
    
    % max inde
    cnt = cnt+1;
    subplottight(r,c,o(cnt),0.15);
    cpg1 = cell(length(compSets),1);
    cpg2 = cell(length(compSets),1);
    for compitr = 1:length(compSets)
        ind = find(mxidx(randIdx) == compitr);
        pg1_ = pg1(randIdx);
        pg2_ = pg2(randIdx);
        cpg1{compitr} = pg1_(ind);
        cpg2{compitr} = pg2_(ind);
    end
    cvec = vecInc(1:length(compSets),cellfun(@numel,cpg1));
    shind = randsample(1:length(cvec),length(cvec));
    cvec = cvec(shind);
    mpg1 = merge(cpg1,1);
    mpg2 = merge(cpg2,1);
    for i = 1:length(cvec)
        plot(mpg2(shind(i)),mpg1(shind(i)),'.','Color',col{cvec(i)},'MarkerSize',msize);
        hold on
    end
    set(gca,'FontSize',fsize)
    xlim([-6,3])
    ylim([-6,7])
    axis square
    ylabel('Gradient 1[uni-trans]')
    xlabel('-Gradient 2[vis-aud]')
    title(sprintf('Best model:%s',p.sbjID{sbjitr}),'FontSize',fsize)
    
    for compitr = 1:length(compSets)
        cnt = cnt+1;
        subplottight(r,c,o(cnt),0.15);
        ind = find(mxidx(randIdx) == compitr);
        pg1_ = pg1(randIdx);
        pg2_ = pg2(randIdx);
        dscatter(pg2_(ind),pg1_(ind),'marker','.','msize',msize);
        set(gca,'FontSize',fsize)
        xlim([-6,3])
        ylim([-6,7])
        axis square
        ylabel('Gradient 1[uni-trans]')
        xlabel('-Gradient 2[vis-aud]')
        title(sprintf('%s:%s',compSets{compitr},p.sbjID{sbjitr}),'FontSize',fsize)
        ffine(h)
    end
    
end

% pooled
% max index
accp = cell(nSbj,1);
pg1p = cell(nSbj,1);
pg2p = cell(nSbj,1);
mxidxp = cell(nSbj,1);
for sbjitr = 1:nSbj
    acc = merge(accs(:,sbjitr),2);
    sig = any(merge(sigs(:,sbjitr),2),2);
    [mx,mxidx] = max(acc,[],2);
    
    % select only signif voxel
    accp{sbjitr} = acc(sig);
    mxidxp{sbjitr} = mxidx(sig);
    pg1p{sbjitr} = pgMaps.(p.sbjID{sbjitr}).map{1}(sig);
    pg2p{sbjitr} = -pgMaps.(p.sbjID{sbjitr}).map{2}(sig);
end
acc = merge(accp,1);
pg1 = merge(pg1p,1);
pg2 = merge(pg2p,1);
mxidx = merge(mxidxp,1);
nSigVox = length(acc);
rand('seed',1);
randIdx = randsample(nSigVox,nPlotVoxel*nSbj);

cpg1 = cell(length(compSets),1);
cpg2 = cell(length(compSets),1);
for compitr = 1:length(compSets)
    ind = find(mxidx(randIdx) == compitr);
    pg1_ = pg1(randIdx);
    pg2_ = pg2(randIdx);
    cpg1{compitr} = pg1_(ind);
    cpg2{compitr} = pg2_(ind);
end
cvec = vecInc(1:length(compSets),cellfun(@numel,cpg1));
shind = randsample(1:length(cvec),length(cvec));
cvec = cvec(shind);
mpg1 = merge(cpg1,1);
mpg2 = merge(cpg2,1);

cnt = cnt+1;
subplottight(r,c,o(cnt),0.15);
for i = 1:length(cvec)
    plot(mpg2(shind(i)),mpg1(shind(i)),'.','Color',col{cvec(i)},'MarkerSize',msize);
    hold on
end
set(gca,'FontSize',fsize)
xlim([-6,3])
ylim([-6,7])
axis square
ylabel('Gradient 1[uni-trans]')
xlabel('-Gradient 2[vis-aud]')
title(sprintf('Best model:Pooled'),'FontSize',fsize)

for compitr = 1:length(compSets)
    cnt = cnt+1;
    subplottight(r,c,o(cnt),0.15);
    ind = find(mxidx(randIdx) == compitr);
    pg1_ = pg1(randIdx);
    pg2_ = pg2(randIdx);
    dscatter(pg2_(ind),pg1_(ind),'marker','.','msize',msize);
    set(gca,'FontSize',fsize)
    xlim([-6,3])
    ylim([-6,7])
    axis square
    ylabel('Gradient 1[uni-trans]')
    xlabel('-Gradient 2[vis-aud]')
    title(sprintf('%s:Pooled',compSets{compitr}),'FontSize',fsize)
    
end

suptitle(sprintf('Distributions of best predicted voxels'));
fname = [p.figdir,'figure5G.pdf'];
pathtitle(fname)
fprintf('%s\n',savprint(h,fname));

%% draw histograms for figure 5G
compSets = {'vision','semantic','category'};


% display settings
fsize = 4;
r = 10;
c = 8;
[r,c,o] = setrc2(r*c,'ltd',[r,c],[1,1]);
thval = 0.111;
nbins = 40;


close all
h = ffigure;
cnt = 0;

accs = cell(length(compSets),nSbj);
sigs = cell(length(compSets),nSbj);
for sbjitr = 1:nSbj
    useIdx = pgMaps.(p.sbjID{sbjitr}).cortexIdx;
    for compitr = 1:length(compSets)
        accs{compitr,sbjitr} = encRes.(compSets{compitr}).pred_acc{sbjitr}(useIdx)';
        sigs{compitr,sbjitr} = accs{compitr,sbjitr} >= thval;
    end
end

for sbjitr = 1:nSbj
    acc = merge(accs(:,sbjitr),2);
    sig = any(merge(sigs(:,sbjitr),2),2);
    [mx,mxidx] = max(acc,[],2);
    
    % select only signif voxel
    mxidx = mxidx(sig);
    pg1 = pgMaps.(p.sbjID{sbjitr}).map{1}(sig);
    pg2 = -pgMaps.(p.sbjID{sbjitr}).map{2}(sig);
    
    
    % max inde
    cpg1 = cell(length(compSets),1);
    cpg2 = cell(length(compSets),1);
    for compitr = 1:length(compSets)
        ind = find(mxidx == compitr);
        cpg1{compitr} = pg1(ind);
        cpg2{compitr} = pg2(ind);
    end
    
    
    % gradient 1
    for compitr = 1:length(compSets)
        cnt = cnt+1;
        subplottight(r,c,o(cnt),0.15);
        bins = linspace(-6,7,nbins);
        [n1,x1] = hist(cpg1{compitr},bins);
        hh = bar(x1,n1,'hist');
        set(hh,'FaceColor',col{compitr},'EdgeColor',[1,1,1])
        set(gca,'FontSize',fsize)
        xlim([-6,7])
        title(sprintf('Grad1:%s:%s',p.sbjID{sbjitr},compSets{compitr}),'FontSize',fsize)
        %ffine(h)
    end
    cnt = cnt+1;
    subplottight(r,c,o(cnt),0.15);
    bins = linspace(-6,7,nbins);
    for compitr = 1:length(compSets)
        val = merge(cpg1((end-compitr+1):-1:1));
        [n1,x1] = hist(val,bins);
        hh = bar(x1,n1,'hist');
        set(hh,'FaceColor',col{end-compitr+1},'EdgeColor',[1,1,1])
        hold on
    end
    set(gca,'FontSize',fsize)
    xlim([-6,7])
    title(sprintf('Grad1:%s:%s',p.sbjID{sbjitr},'All'),'FontSize',fsize)
    %ffine(h)
    
    % gradient 2
    for compitr = 1:length(compSets)
        cnt = cnt+1;
        subplottight(r,c,o(cnt),0.15);
        bins = linspace(-6,3,nbins);
        [n1,x1] = hist(cpg2{compitr},bins);
        hh = bar(x1,n1,'hist');
        set(hh,'FaceColor',col{compitr},'EdgeColor',[1,1,1])
        set(gca,'FontSize',fsize)
        xlim([-6,3])
        title(sprintf('Grad2:%s:%s',p.sbjID{sbjitr},compSets{compitr}),'FontSize',fsize)
        %ffine(h)
    end
    cnt = cnt+1;
    subplottight(r,c,o(cnt),0.15);
    bins = linspace(-6,3,nbins);
    for compitr = 1:length(compSets)
        val = merge(cpg2((end-compitr+1):-1:1));
        [n1,x1] = hist(val,bins);
        hh = bar(x1,n1,'hist');
        set(hh,'FaceColor',col{end-compitr+1},'EdgeColor',[1,1,1])
        hold on
    end
    set(gca,'FontSize',fsize)
    xlim([-6,3])
    title(sprintf('Grad2:%s:%s',p.sbjID{sbjitr},'All'),'FontSize',fsize)
    %ffine(h)
end

% pooled
% max index
accp = cell(nSbj,1);
pg1p = cell(nSbj,1);
pg2p = cell(nSbj,1);
mxidxp = cell(nSbj,1);
for sbjitr = 1:nSbj
    acc = merge(accs(:,sbjitr),2);
    sig = any(merge(sigs(:,sbjitr),2),2);
    [mx,mxidx] = max(acc,[],2);
    
    % select only signif voxel
    accp{sbjitr} = acc(sig);
    mxidxp{sbjitr} = mxidx(sig);
    pg1 = pgMaps.(p.sbjID{sbjitr}).map{1};
    pg2 = -pgMaps.(p.sbjID{sbjitr}).map{2};
    pg1p{sbjitr} = pg1(sig);
    pg2p{sbjitr} = pg2(sig);
    
end
acc = merge(accp,1);
pg1 = merge(pg1p,1);
pg2 = merge(pg2p,1);
mxidx = merge(mxidxp,1);

cpg1 = cell(length(compSets),1);
cpg2 = cell(length(compSets),1);
for compitr = 1:length(compSets)
    ind = find(mxidx == compitr);
    cpg1{compitr} = pg1(ind);
    cpg2{compitr} = pg2(ind);
end


% gradient 1
for compitr = 1:length(compSets)
    cnt = cnt+1;
    subplottight(r,c,o(cnt),0.15);
    bins = linspace(-6,7,nbins);
    [n1,x1] = hist(cpg1{compitr},bins);
    hh = bar(x1,n1,'hist');
    set(hh,'FaceColor',col{compitr},'EdgeColor',[1,1,1])
    set(gca,'FontSize',fsize)
    xlim([-6,7])
    title(sprintf('Grad1:%s:%s:','Pooled',compSets{compitr}),'FontSize',fsize)
    %ffine(h)
end
cnt = cnt+1;
subplottight(r,c,o(cnt),0.15);
bins = linspace(-6,7,nbins);
for compitr = 1:length(compSets)
    val = merge(cpg1((end-compitr+1):-1:1));
    [n1,x1] = hist(val,bins);
    hh = bar(x1,n1,'hist');
    set(hh,'FaceColor',col{end-compitr+1},'EdgeColor',[1,1,1])
    hold on
end
set(gca,'FontSize',fsize)
xlim([-6,7])
title(sprintf('Grad1:%s:%s','Pooled','All'),'FontSize',fsize)
%ffine(h)

% gradient 2
for compitr = 1:length(compSets)
    cnt = cnt+1;
    subplottight(r,c,o(cnt),0.15);
    bins = linspace(-6,3,nbins);
    [n1,x1] = hist(cpg2{compitr},bins);
    hh = bar(x1,n1,'hist');
    set(hh,'FaceColor',col{compitr},'EdgeColor',[1,1,1])
    set(gca,'FontSize',fsize)
    xlim([-6,3])
    title(sprintf('Grad2:%s:%s:','Pooled',compSets{compitr}),'FontSize',fsize)
    %ffine(h)
end
cnt = cnt+1;
subplottight(r,c,o(cnt),0.15);
bins = linspace(-6,3,nbins);
for compitr = 1:length(compSets)
    val = merge(cpg2((end-compitr+1):-1:1));
    [n1,x1] = hist(val,bins);
    hh = bar(x1,n1,'hist');
    set(hh,'FaceColor',col{end-compitr+1},'EdgeColor',[1,1,1])
    hold on
end
set(gca,'FontSize',fsize)
xlim([-6,3])
title(sprintf('Grad2:%s:%s','Pooled','All'),'FontSize',fsize)
%ffine(h)


suptitle(sprintf('Histograms of best predicted voxels'));
fname = [p.figdir,'figure5G_hist.pdf'];
pathtitle(fname)
fprintf('%s\n',savprint(h,fname));


%% draw figure 5H
compSets = {'semantic','category'};

axrange = [1,5,10];
nbins = 10;
range = [0,0.15];

% display settings
r = 4;
c = 9;
[r,c,o] = setrc2(r*c,'ltr',[r,c],[1,2]);
fsize = 6;
msize = 3;%1.5;

close all
h = ffigure;
cnt = 0;


accs = cell(length(compSets),nSbj);
for sbjitr = 1:nSbj
    useIdx = pgMaps.(p.sbjID{sbjitr}).cortexIdx;
    for compitr = 1:length(compSets)
        accs{compitr,sbjitr} = encRes.(compSets{compitr}).pred_acc{sbjitr}(useIdx)';
    end
end

BinIdx = cell(nSbj,2);
for sbjitr = 1:nSbj
    
    pg1 = pgMaps.(p.sbjID{sbjitr}).map{1};
    pg2 = -pgMaps.(p.sbjID{sbjitr}).map{2};
    
    % select equivalent numbers of voxels for each bin
    [sorted,ord] = sort(pg1);
    BinIdx{sbjitr,1} = sorted(floor(linspace(1,length(pg1),nbins+1)))';
    [sorted,ord] = sort(pg2);
    BinIdx{sbjitr,2} = sorted(floor(linspace(1,length(pg2),nbins+1)))';
end

% category semantic difference
pMat = zeros(nSbj,nbins);
catsemdiffMat = cell(2,1);
acc1Mat_mu = cell(nSbj,1);
acc1Mat_ci = cell(nSbj,1);
acc2Mat_mu = cell(nSbj,1);
acc2Mat_ci = cell(nSbj,1);

for sbjitr = 1:nSbj
    pg1 = pgMaps.(p.sbjID{sbjitr}).map{1};
    useIdx = pgMaps.(p.sbjID{sbjitr}).cortexIdx;
    
    acc1 = encRes.(compSets{ismember(compSets,compSets{1})}).pred_acc{sbjitr}(useIdx)';
    acc2 = encRes.(compSets{ismember(compSets,compSets{2})}).pred_acc{sbjitr}(useIdx)';
    
    
    for pgitr = 1:nbins
        % pg1
        idx = (pg1 > BinIdx{sbjitr,1}(pgitr)) & pg1 <= BinIdx{sbjitr,1}(pgitr+1);
        [acc1Mat_ci{sbjitr,1}(pgitr),acc1Mat_mu{sbjitr,1}(pgitr)] = ciestim3(acc1(idx));
        [acc2Mat_ci{sbjitr,1}(pgitr),acc2Mat_mu{sbjitr,1}(pgitr)] = ciestim3(acc2(idx));
        
        [sig,pMat(sbjitr,pgitr)] = ttest(r2z(acc1(idx)),r2z(acc2(idx)));
        catsemdiffMat{1}(sbjitr,pgitr) = sum(acc1(idx)-acc2(idx));
        catsemdiffMat{2}(sbjitr,pgitr) = sum(acc2(idx)-acc1(idx));
        
    end
end
% draw plot
mcompcorr = nSbj*nbins;
th = 0.001;
th_act = th/mcompcorr;
for sbjitr = 1:nSbj
    cnt = cnt+1;
    subplottight(r,c,o(cnt),0.18);
    for ix = 1:length(compSets)
        switch compSets{ix}
            case 'semantic'
                accMat_mu = acc1Mat_mu;
            case 'category'
                accMat_mu = acc2Mat_mu;
        end
        acc = accMat_mu{sbjitr};
        
        seq = 1:length(acc);
        hh = plot(seq,acc,'-');
        set(hh,'Color',plotColor(ix,length(compSets),'rbemo'))
        hold on
        
        sig = pMat(sbjitr,:) < th_act &  catsemdiffMat{ix}(sbjitr,:) > 0;
        nsac = accMat_mu{sbjitr}(~sig)';
        hh = plot(seq(~sig),nsac,'o');
        set(hh,'Color',plotColor(ix,length(compSets),'rbemo'),'MarkerFaceColor',plotColor(ix,length(compSets),'rbemo'),'MarkerSize',msize)
        AnnotationOff(hh,1:length(hh));
        
        sac = accMat_mu{sbjitr}(sig)';
        hh = plot(seq(sig),sac,'o');
        set(hh,'Color',plotColor(ix,length(compSets),'rbemo'),'MarkerFaceColor',[1,1,1],'MarkerSize',msize)
        AnnotationOff(hh,1:length(hh));
    end
    
    set(gca,'FontSize',fsize)
    ylabel(sprintf('Prediction performance (r)'))
    axname(axrange,1,axrange)
    xlabel('PG level')
    ylim(range)
    xlim([0,nbins+1])
    hline(0,'-k')
    legend(compSets,'Location','SouthWest')
    ffine(h)
    title(sprintf('%s',p.sbjID{sbjitr}),'FontSize',fsize)
end


suptitle(sprintf('Mean encoding accuracy along first PG axis'));
fname = [p.figdir,'figure5H.pdf'];
pathtitle(fname)
fprintf('%s\n',savprint(h,fname));

%%
close all
