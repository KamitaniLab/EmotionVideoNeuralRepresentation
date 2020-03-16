function draw_figure3ABC_S1C(p,decRes)
%
% This code is for drawing figure 3ABC and S1C
%
%
%% settings
compPairs = {'category','dimension'};

% roi parameters
hcpIdx = ~cellfun(@isempty,strfind(p.roiDescrip,'hcp180'));
subIdx = ismember(p.roiDescrip,{'thalamus','hyppocampus','hypothalamus','pallidum','brainstem','caudate','putamen','nuc_accum','amygdala','cerebellum'});
roiIdx = hcpIdx|subIdx;
nhcp = sum(hcpIdx);
nsub = sum(subIdx);
nRois = length(p.roiDescrip);


% get color index for hcp360 rois
tmp = load([p.rootPath,'code/util/hcp360colors.mat']);
hcpcols = tmp.vrank;

colcat = [216,31,73]/255;
coldim = [250,120,210]/255;

nSbj = length(p.sbjID);

%% draw figure 3A and S1C

% display settings
r = 5;
c = 5;
[r,c,o] = setrc2(r*c,'ltr',[r,c],[1]);
fsize = 8;
msize = 4;%1.5;
msize_sub = 1.2;

close all
h = ffigure;
cmap3('ibg4');
cnt = 0;


Mu = cell(nSbj,1);
for sbjitr = 1:nSbj
    % summary
    mu = zeros(nRois,length(compPairs));
    for ix = 1:length(compPairs)
        md = decRes.(compPairs{ix}){sbjitr}.mRoiDec;
        mu(:,ix) = mean(md.iden_acc_all,1);
    end
    Mu{sbjitr} = mu;
    
    cnt = cnt+1;
    subplottight(r,c,o(cnt),0.15);
    mus = mu(roiIdx,:);
    [s,ord] = sort(mean(mus,2),'ascend');
    mu_sorted = mus(ord,:);
    
    % draw significant ROIs
    for ix = 1:length(mu_sorted)
        hold on
        % cortex
        if ord(ix) <= nhcp
            plot(mu_sorted(ix,1),mu_sorted(ix,2),'.','MarkerEdgeColor',plotColor(hcpcols(ord(ix)),(nhcp+nsub)/2,'Spectral_r'),'MarkerSize',msize,'LineWidth',0.1);
        elseif (ord(ix) == nhcp+8)||(ord(ix) == nhcp+9) % +x
            % subcortex
            plot(mu_sorted(ix,1),mu_sorted(ix,2),plotMarker2(ord(ix)-nhcp),...
                'MarkerEdgeColor',[1,1,1]*0.2,...
                'MarkerFaceColor',[1,1,1]*0.2,...
                'MarkerSize',1.6,'LineWidth',0.2);
        else
            % subcortex
            plot(mu_sorted(ix,1),mu_sorted(ix,2),plotMarker2(ord(ix)-nhcp),...
                'MarkerEdgeColor',[1,1,1]*0.2,...
                'MarkerFaceColor',[1,1,1]*0.2,...
                'MarkerSize',msize_sub,'LineWidth',0.1);
        end
    end
    hold on
    axis square
    xlabel(compPairs{1})
    ylabel(compPairs{2})
    xlim([40,100])
    ylim([40,100])
    odline(axis)
    xlim([45,80])
    ylim([45,80])
    hline(50,'--k');
    vline(50,'--k');
    ffine(h)
    title(sprintf('%s vs. %s: %s',compPairs{1},compPairs{2},p.sbjID{sbjitr}),'FontSize',fsize)
    
end

% mean
cnt = cnt+1;
subplottight(r,c,o(cnt),0.15);
mu = cellmean(Mu);
mus = mu(roiIdx,:);
[s,ord] = sort(mean(mus,2),'ascend');
mu_sorted = mus(ord,:);


for ix = 1:length(mu_sorted)
    hold on
    % cortex
    if ord(ix) <= nhcp
        plot(mu_sorted(ix,1),mu_sorted(ix,2),'.','MarkerEdgeColor',plotColor(hcpcols(ord(ix)),(nhcp+nsub)/2,'Spectral_r'),'MarkerSize',msize,'LineWidth',0.1);
    elseif (ord(ix) == nhcp+8)||(ord(ix) == nhcp+9) % +x
        % subcortex
        plot(mu_sorted(ix,1),mu_sorted(ix,2),plotMarker2(ord(ix)-nhcp),...
            'MarkerEdgeColor',[1,1,1]*0.2,...
            'MarkerFaceColor',[1,1,1]*0.2,...
            'MarkerSize',1.6,'LineWidth',0.2);
    else
        % subcortex
        plot(mu_sorted(ix,1),mu_sorted(ix,2),plotMarker2(ord(ix)-nhcp),...
            'MarkerEdgeColor',[1,1,1]*0.2,...
            'MarkerFaceColor',[1,1,1]*0.2,...
            'MarkerSize',msize_sub,'LineWidth',0.1);
    end
end
hold on

axis square
xlabel(compPairs{1})
ylabel(compPairs{2})
xlim([40,100])
ylim([40,100])
odline(axis)
xlim([45,80])
ylim([45,80])
hline(50,'--k');
vline(50,'--k');
ffine(h)
title(sprintf('%s vs. %s: %s',compPairs{1},compPairs{2},'Average'),'FontSize',fsize)


suptitle(sprintf('Identification accuracy via region-wise decoding analyses'));
fname = [p.figdir,'figure3A_S1C.pdf'];
pathtitle(fname)
fprintf('%s\n',savprint(h,fname));


%% draw figures 3B

% summarize idnetification accuracy
iden_cr_mu = cell(length(compPairs),nSbj);
iden_cr_ci = cell(length(compPairs),nSbj);
for sbjitr = 1:nSbj
    for scoritr = 1:length(compPairs)
        fset1 = compPairs{scoritr};
        ed = decRes.(fset1){sbjitr}.ensDec;
        acc1 = ed.idenacc_mean;
        
        [ci,mu] = ciestim3(acc1,1,0.99,'onesided');
        iden_cr_mu{scoritr,sbjitr} = mu;
        iden_cr_ci{scoritr,sbjitr} = ci;
        
    end
end


% display settings
r = 5;
c = 5;
[r,c,o] = setrc2(r*c,'ltr',[r,c],[2,1]);
fsize = 8;
msize = 2;%1.5;

close all
h = ffigure;
cmap3('ibg4');
cnt = 0;


cnt = cnt + 1;
subplottight(r,c,o(cnt),0.1);
% summary
mu = zeros(nSbj,length(compPairs));
ci = zeros(nSbj,length(compPairs));
for ix = 1:length(compPairs)
    for sbjitr = 1:nSbj
        scoreidx = ismember(compPairs,compPairs{ix});
        mu(sbjitr,ix) = iden_cr_mu{scoreidx,sbjitr}*100;
        ci(sbjitr,ix) = iden_cr_ci{scoreidx,sbjitr}*100;
    end
end

hold on
bar(mean(mu,1)','edgecolor','none')
errorbar_h(mean(mu,1)',ciestim3(mu,1,0.99,'onesided')','.k');
dif = 0.2;
xloc = repmat((1:length(compPairs))+dif,size(mu,1),1);
xloc = xloc+repmat(randn(size(xloc,1),1)/50,1,2);
plot(xloc',mu','-o','Color',[1,1,1]*0.4,'MarkerEdgeColor',[1,1,1]*0.2,'MarkerFaceColor',[1,1,1]*0.4,'MarkerSize',msize);
set(gca,'FontSize',fsize)
ylabel(sprintf('Identification accuracy (%%)'))
axname(compPairs,1)
ylim([40,100])
xlim([0,6])
xticklabel_rotate([],45)
hline(50,'-k')
hh = hline(40:10:100,'-k');
set(hh,'Color',[1,1,1]*0.8)
axis square
ffine(h)
suptitle(sprintf('Mean identification accuracy via decoding analyses'));
fname = [p.figdir,'figure3B.pdf'];
pathtitle(fname)
fprintf('%s\n',savprint(h,fname));


%% draw figure 3C

% display settings
c = 5;
r = 6;
[r,c,o] = setrc2(r*c,'ltd',[r,c],[1]);
fsize = 8;
sbjs = {p.sbjID{:},'Averaged'};

AccROIs_cat = cell(length(p.sbjID),1);
AccROIs_dim = cell(length(p.sbjID),1);
for sbjitr = 1:length(p.sbjID)
    catRes = decRes.category{sbjitr};
    dimRes = decRes.dimension{sbjitr};
    
    AccROIs_dim{sbjitr} = dimRes.ensDec.idenacc_mean;
    AccROIs_cat{sbjitr} = catRes.ensDec.idenacc_mean;
    
end
clear Cnt
h = ffigure;
cnt = 0;

for sbjitr = 1:(nSbj+1)
    if sbjitr > nSbj
        acc = cellmean(AccROIs_cat)*100;
        acd = cellmean(AccROIs_dim)*100;
        rang = [0,250];
    else
        acc = AccROIs_cat{sbjitr}*100;
        acd = AccROIs_dim{sbjitr}*100;
        rang = [0,550];
    end
    
    % category
    [n0,x0] = hist([acc;acd],40);
    
    cnt = cnt+1;
    subplottight(r,c,o(cnt),0.1);
    [n1,x1] = hist(acc,x0);
    hh = bar(x1,n1,'hist');
    set(gca,'FontSize',fsize)
    set(hh,'facecolor',colcat,'edgecolor',[1,1,1])
    xlim([0,100])
    ylim(rang)
    vline(50,'-k')
    vline(mean(acc),'--k')
    axname(0:10:100,1,0:10:100)
    xlabel('Identification accuracy (%)')
    ylabel('Frequency')
    title(sprintf('Category:%s',sbjs{sbjitr}))
    ffine(h)
    
    % dimension
    cnt = cnt+1;
    subplottight(r,c,o(cnt),0.1);
    [n1,x1] = hist(acd,x0);
    hh = bar(x1,n1,'hist');
    set(gca,'FontSize',fsize)
    set(hh,'facecolor',coldim,'edgecolor',[1,1,1])
    xlim([0,100])
    ylim(rang)
    vline(50,'-k')
    vline(mean(acd),'--k')
    axname(0:10:100,1,0:10:100)
    xlabel('Identification accuracy (%)')
    ylabel('Frequency')
    title(sprintf('Dimension:%s',sbjs{sbjitr}))
    ffine(h)
    
end
suptitle(sprintf('Histogram of decoding accuracy for individual videos'));
fname = [p.figdir,'figure3C.pdf'];
pathtitle(fname)
fprintf('%s\n',savprint(h,fname));

%%
close all
