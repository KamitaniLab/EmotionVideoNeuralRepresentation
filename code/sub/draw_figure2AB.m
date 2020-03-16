function draw_figure2AB(p,decRes)
%
% This code is for drawing figure 2AB
%
%
%% settings
% score index
dimIdx = 1:14;
catIdx = 15:48;

% get roi Idx
hcpIdx = ~cellfun(@isempty,strfind(p.roiDescrip,'hcp180'));
subIdx = ismember(p.roiDescrip,{'thalamus','hyppocampus','hypothalamus','pallidum','brainstem','caudate','putamen','nuc_accum','amygdala','cerebellum'});

roiIdx = hcpIdx|subIdx;
nRois = sum(roiIdx);
nhcp = sum(hcpIdx);

% get color index for hcp360 rois
tmp = load([p.rootPath,'code/util/hcp360colors.mat']);
hcpcols = tmp.vrank;

colcat = [216,31,73]/255;
coldim = [250,120,210]/255;

nSbj = length(p.sbjID);

% get feature names
load([p.featdir,'category.mat']);
catname = L.featname';
load([p.featdir,'dimension.mat']);
dimname = L.featname';

%% get figure plot parameters

AccEns = cell(nSbj,1);
AccROIs = cell(nSbj,1);
AccROIsCV = cell(nSbj,1);
for sbjitr = 1:nSbj
    catRes = decRes.category{sbjitr};
    dimRes = decRes.dimension{sbjitr};
    
    AccEns{sbjitr} = [dimRes.ensDec.profile_acc;catRes.ensDec.profile_acc];
    AccROIs{sbjitr} = [dimRes.mRoiDec.profile_acc_all(:,roiIdx);catRes.mRoiDec.profile_acc_all(:,roiIdx)];
    AccROIsCV{sbjitr} = [dimRes.mRoiDec.profile_cvacc_all(:,roiIdx);catRes.mRoiDec.profile_cvacc_all(:,roiIdx)];
end

% get plot params ==
clear SIG
xd = cell(nSbj+1,1);
yd = cell(nSbj+1,1);
o_dim = cell(nSbj+1,1);
o_cat = cell(nSbj+1,1);
o_both = cell(nSbj+1,1);
for sbjitr = 1:nSbj
    [s_dim,o_dim{sbjitr}] = sort(AccEns{sbjitr}(dimIdx,:),'descend');
    [s_cat,o_cat{sbjitr}] = sort(AccEns{sbjitr}(catIdx,:),'descend');
    o_catd = o_cat{sbjitr}+length(o_dim{sbjitr});
    o_both{sbjitr} = [o_catd;o_dim{sbjitr}];
    
    hh = sinaplot(AccROIs{sbjitr}(o_both{sbjitr},:)');
    
    xd{sbjitr} = get(hh,'XData');
    yd{sbjitr} = get(hh,'YData');
end

% Average
cAccROIs = cellmean(AccEns);
[s_dim,o_dim{sbjitr+1}] = sort(cAccROIs(dimIdx,:),'descend');
[s_cat,o_cat{sbjitr+1}] = sort(cAccROIs(catIdx,:),'descend');
o_catd = o_cat{sbjitr+1}+length(o_dim{sbjitr+1});
o_both{sbjitr+1} = [o_catd;o_dim{sbjitr+1}];
cAccROIs = cellmean(AccROIs);

hh = sinaplot(cAccROIs(o_both{sbjitr+1},:)');
xd{sbjitr+1} = get(hh,'XData');
yd{sbjitr+1} = get(hh,'YData');
close all
% end get figure params ==


%% draw figures

% display settings
c = 2;
r = 6;
[r,c,o] = setrc2(r*c,'ltr',[r,c],[1,0]);
fsize = 6;
msize = 2;
msize_sub = 0.6;
msize_nsig = 0.4;
thval = 0.095;
cnt = 0;

close all
h = ffigure;
SIG = cell(nSbj,1);
for sbjitr = 1:nSbj
    cnt = cnt+1;
    subplottight(r,c,o(cnt),0.1);
    
    for i = 1:length(xd{sbjitr})
        sig1 = AccROIsCV{sbjitr}(o_both{sbjitr}(i),:)' > thval;
        hold on
        if i > length(catIdx)
            sigIdx = find(sig1);
            nonsigIdx = find(~sig1);
            [s,ord] = sort(yd{sbjitr}{i}(sigIdx),'ascend');
            sig1Idx_sorted = sigIdx(ord);
            hcpcolssig = hcpcols(sig1(1:nhcp));
            % draw non significant ROIs
            % cortex
            hold on
            plot(xd{sbjitr}{i}(nonsigIdx(nonsigIdx<=nhcp)),yd{sbjitr}{i}(nonsigIdx(nonsigIdx<=nhcp)),'o','MarkerEdgeColor',[1,1,1]*0.8,'MarkerSize',msize_nsig,'LineWidth',0.1);
            % subcortex
            if any(nonsigIdx(nonsigIdx>nhcp))
                subnonsigIdx = nonsigIdx(nonsigIdx>nhcp);
                for ixxx = 1:length(subnonsigIdx)
                    plot(xd{sbjitr}{i}(subnonsigIdx(ixxx)),yd{sbjitr}{i}(subnonsigIdx(ixxx)),plotMarker2(subnonsigIdx(ixxx)-nhcp),'MarkerEdgeColor',[1,1,1]*0.8,'MarkerSize',msize_nsig,'LineWidth',0.1);
                end
            end
            
            % draw significant ROIs
            for ix = 1:length(sig1Idx_sorted)
                % cortex
                if sig1Idx_sorted(ix) <= nhcp
                    plot(xd{sbjitr}{i}(sig1Idx_sorted(ix)),yd{sbjitr}{i}(sig1Idx_sorted(ix)),'.','MarkerEdgeColor',plotColor(hcpcolssig(ord(ix)),nRois/2,'Spectral_r'),'MarkerSize',msize,'LineWidth',0.1);
                elseif (sig1Idx_sorted(ix) == nhcp+8) || (sig1Idx_sorted(ix) == nhcp+9) % +x
                    % subcortex
                    plot(xd{sbjitr}{i}(sig1Idx_sorted(ix)),yd{sbjitr}{i}(sig1Idx_sorted(ix)),plotMarker2(sig1Idx_sorted(ix)-nhcp),...
                        'MarkerEdgeColor',[1,1,1]*0.2,...
                        'MarkerFaceColor',[1,1,1]*0.2,...
                        'MarkerSize',0.8,'LineWidth',0.2);
                else
                    % subcortex
                    plot(xd{sbjitr}{i}(sig1Idx_sorted(ix)),yd{sbjitr}{i}(sig1Idx_sorted(ix)),plotMarker2(sig1Idx_sorted(ix)-nhcp),...
                        'MarkerEdgeColor',[1,1,1]*0.2,...
                        'MarkerFaceColor',[1,1,1]*0.2,...
                        'MarkerSize',msize_sub,'LineWidth',0.1);
                end
            end
            
        else % = category
            sigIdx = find(sig1);
            nonsigIdx = find(~sig1);
            [s,ord] = sort(yd{sbjitr}{i}(sigIdx),'ascend');
            sig1Idx_sorted = sigIdx(ord);
            hcpcolssig = hcpcols(sig1(1:nhcp));
            % draw non significant ROIs
            % cortex
            hold on
            plot(xd{sbjitr}{i}(nonsigIdx(nonsigIdx<=nhcp)),yd{sbjitr}{i}(nonsigIdx(nonsigIdx<=nhcp)),'o','MarkerEdgeColor',[1,1,1]*0.8,'MarkerSize',msize_nsig,'LineWidth',0.1);
            % subcortex
            if any(nonsigIdx(nonsigIdx>nhcp))
                subnonsigIdx = nonsigIdx(nonsigIdx>nhcp);
                for ixxx = 1:length(subnonsigIdx)
                    plot(xd{sbjitr}{i}(subnonsigIdx(ixxx)),yd{sbjitr}{i}(subnonsigIdx(ixxx)),plotMarker2(subnonsigIdx(ixxx)-nhcp),'MarkerEdgeColor',[1,1,1]*0.8,'MarkerSize',msize_nsig,'LineWidth',0.1);
                end
            end
            % draw significant ROIs
            for ix = 1:length(sig1Idx_sorted)
                % cortex
                if sig1Idx_sorted(ix) <= nhcp
                    plot(xd{sbjitr}{i}(sig1Idx_sorted(ix)),yd{sbjitr}{i}(sig1Idx_sorted(ix)),'.','MarkerEdgeColor',plotColor(hcpcolssig(ord(ix)),nRois/2,'Spectral_r'),'MarkerSize',msize,'LineWidth',0.1);
                elseif (sig1Idx_sorted(ix) == nhcp+8) || (sig1Idx_sorted(ix) == nhcp+9) % +x
                    % subcortex
                    plot(xd{sbjitr}{i}(sig1Idx_sorted(ix)),yd{sbjitr}{i}(sig1Idx_sorted(ix)),plotMarker2(sig1Idx_sorted(ix)-nhcp),...
                        'MarkerEdgeColor',[1,1,1]*0.2,...
                        'MarkerFaceColor',[1,1,1]*0.2,...
                        'MarkerSize',0.8,'LineWidth',0.2);
                else
                    % subcortex
                    plot(xd{sbjitr}{i}(sig1Idx_sorted(ix)),yd{sbjitr}{i}(sig1Idx_sorted(ix)),plotMarker2(sig1Idx_sorted(ix)-nhcp),...
                        'MarkerEdgeColor',[1,1,1]*0.2,...
                        'MarkerFaceColor',[1,1,1]*0.2,...
                        'MarkerSize',msize_sub,'LineWidth',0.1);
                end
            end
            
        end
        SIG{sbjitr}(:,i) = AccROIsCV{sbjitr}(i,:)' > thval;
    end
    
    
    hold on
    mm = AccEns{sbjitr}(o_both{sbjitr},:);
    mm(isnan(mm)) = 0;
    
    plot(1:length(catIdx),mm(1:length(catIdx)),'s','MarkerEdgeColor',colcat,'MarkerSize',msize)
    plot((1:length(dimIdx))+length(catIdx),mm(length(catIdx)+1:end),'s','MarkerEdgeColor',coldim,'MarkerSize',msize)
    ylim([-0.1,0.6])
    set(gca,'FontSize',fsize)
    hh = hline(0.1:0.1:0.7,'-k');
    set(hh,'Color',[1,1,1]*0.8)
    hline(0,'-k');
    vline(length(catIdx)+0.5,'-k');
    ylabel('Corr.Coeff.')
    axname(strrep([catname(o_cat{sbjitr}),dimname(o_dim{sbjitr})],'_',' '),1)
    xticklabel_rotate([],45)
    title(sprintf('%s',p.sbjID{sbjitr}))
    ffine(h)
end



% draw figures for averaged
sigsum = cellsum(SIG)';
cnt = cnt+1;
subplottight(r,c,o(cnt),0.1);
hold on

% ensemble accuracy
[ci,mu] = cellci(AccEns);
mu = mu(o_both{sbjitr+1});
ci = ci(o_both{sbjitr+1});
for i = 1:length(mu)
    if i > length(catIdx)
        hh = errorbar_h(i,mu(i),ci(i),'.');
        col = coldim;
        set(hh,'Color',col,'MarkerEdgeColor',col,'MarkerFaceColor',col,'LineWidth',0.5)
        hh = plot(i,mu(i),'s');
        set(hh,'Color',col,'MarkerEdgeColor',col,'LineWidth',1)
    else % = category
        hh = errorbar_h(i,mu(i),ci(i),'.');
        col = colcat;
        set(hh,'Color',col,'MarkerEdgeColor',col,'MarkerFaceColor',col,'LineWidth',0.5)
        hh = plot(i,mu(i),'s');
        set(hh,'Color',col,'MarkerEdgeColor',col,'LineWidth',1)
    end
    set(hh,'MarkerSize',1)
end
ylim([-0.1,0.6])
set(gca,'FontSize',fsize)

for i = 1:length(xd{sbjitr+1})
    hold on
    sig1 = sigsum(o_both{sbjitr+1}(i),:)==nSbj;
    if i > length(catIdx)
        sigIdx = find(sig1);
        nonsigIdx = find(~sig1);
        [s,ord] = sort(yd{sbjitr+1}{i}(sigIdx),'ascend');
        sig1Idx_sorted = sigIdx(ord);
        hcpcolssig = hcpcols(sig1(1:nhcp));
        % draw non significant ROIs
        % cortex
        hold on
        plot(xd{sbjitr+1}{i}(nonsigIdx(nonsigIdx<=nhcp)),yd{sbjitr+1}{i}(nonsigIdx(nonsigIdx<=nhcp)),'o','MarkerEdgeColor',[1,1,1]*0.8,'MarkerSize',msize_nsig,'LineWidth',0.1);
        % subcortex
        if any(nonsigIdx(nonsigIdx>nhcp))
            subnonsigIdx = nonsigIdx(nonsigIdx>nhcp);
            for ixxx = 1:length(subnonsigIdx)
                plot(xd{sbjitr+1}{i}(subnonsigIdx(ixxx)),yd{sbjitr+1}{i}(subnonsigIdx(ixxx)),plotMarker2(subnonsigIdx(ixxx)-nhcp),'MarkerEdgeColor',[1,1,1]*0.8,'MarkerSize',msize_nsig,'LineWidth',0.1);
            end
        end
        % draw significant ROIs
        for ix = 1:length(sig1Idx_sorted)
            % cortex
            if sig1Idx_sorted(ix) <= nhcp
                plot(xd{sbjitr+1}{i}(sig1Idx_sorted(ix)),yd{sbjitr+1}{i}(sig1Idx_sorted(ix)),'.','MarkerEdgeColor',plotColor(hcpcolssig(ord(ix)),nRois/2,'Spectral_r'),'MarkerSize',msize,'LineWidth',0.1);
            elseif (sig1Idx_sorted(ix) == nhcp+8) || (sig1Idx_sorted(ix) == nhcp+9) % +x
                % subcortex
                plot(xd{sbjitr+1}{i}(sig1Idx_sorted(ix)),yd{sbjitr+1}{i}(sig1Idx_sorted(ix)),plotMarker2(sig1Idx_sorted(ix)-nhcp),...
                    'MarkerEdgeColor',[1,1,1]*0.2,...
                    'MarkerFaceColor',[1,1,1]*0.2,...
                    'MarkerSize',0.8,'LineWidth',0.2);
            else
                % subcortex
                plot(xd{sbjitr+1}{i}(sig1Idx_sorted(ix)),yd{sbjitr+1}{i}(sig1Idx_sorted(ix)),plotMarker2(sig1Idx_sorted(ix)-nhcp),...
                    'MarkerEdgeColor',[1,1,1]*0.2,...
                    'MarkerFaceColor',[1,1,1]*0.2,...
                    'MarkerSize',msize_sub,'LineWidth',0.1);
            end
        end
        
    else % = category
        sigIdx = find(sig1);
        nonsigIdx = find(~sig1);
        [s,ord] = sort(yd{sbjitr+1}{i}(sigIdx),'ascend');
        sig1Idx_sorted = sigIdx(ord);
        hcpcolssig = hcpcols(sig1(1:nhcp));
        % draw non significant ROIs
        % cortex
        hold on
        plot(xd{sbjitr+1}{i}(nonsigIdx(nonsigIdx<=nhcp)),yd{sbjitr+1}{i}(nonsigIdx(nonsigIdx<=nhcp)),'o','MarkerEdgeColor',[1,1,1]*0.8,'MarkerSize',msize_nsig,'LineWidth',0.1);
        % subcortex
        if any(nonsigIdx(nonsigIdx>nhcp))
            subnonsigIdx = nonsigIdx(nonsigIdx>nhcp);
            for ixxx = 1:length(subnonsigIdx)
                plot(xd{sbjitr+1}{i}(subnonsigIdx(ixxx)),yd{sbjitr+1}{i}(subnonsigIdx(ixxx)),plotMarker2(subnonsigIdx(ixxx)-nhcp),'MarkerEdgeColor',[1,1,1]*0.8,'MarkerSize',msize_nsig,'LineWidth',0.1);
            end
        end
        % draw significant ROIs
        for ix = 1:length(sig1Idx_sorted)
            % cortex
            if sig1Idx_sorted(ix) <= nhcp
                plot(xd{sbjitr+1}{i}(sig1Idx_sorted(ix)),yd{sbjitr+1}{i}(sig1Idx_sorted(ix)),'.','MarkerEdgeColor',plotColor(hcpcolssig(ord(ix)),nRois/2,'Spectral_r'),'MarkerSize',msize,'LineWidth',0.1);
            elseif (sig1Idx_sorted(ix) == nhcp+8) || (sig1Idx_sorted(ix) == nhcp+9) % +x
                % subcortex
                plot(xd{sbjitr+1}{i}(sig1Idx_sorted(ix)),yd{sbjitr+1}{i}(sig1Idx_sorted(ix)),plotMarker2(sig1Idx_sorted(ix)-nhcp),...
                    'MarkerEdgeColor',[1,1,1]*0.2,...
                    'MarkerFaceColor',[1,1,1]*0.2,...
                    'MarkerSize',0.8,'LineWidth',0.2);
            else
                % subcortex
                plot(xd{sbjitr+1}{i}(sig1Idx_sorted(ix)),yd{sbjitr+1}{i}(sig1Idx_sorted(ix)),plotMarker2(sig1Idx_sorted(ix)-nhcp),...
                    'MarkerEdgeColor',[1,1,1]*0.2,...
                    'MarkerFaceColor',[1,1,1]*0.2,...
                    'MarkerSize',msize_sub,'LineWidth',0.1);
            end
        end
        
    end
end
hh = hline(0.1:0.1:0.7,'-k');
set(hh,'Color',[1,1,1]*0.8)
hline(0,'-k');
vline(length(catIdx)+0.5,'-k');
ylabel('Corr.Coeff.')
axname(strrep([catname(o_cat{sbjitr+1}),dimname(o_dim{sbjitr+1})],'_',' '),1)
xticklabel_rotate([],45)
title(sprintf('%s','Average'))
ffine(h)

suptitle(sprintf('Decoding accuracy for emotion scores'));
fname = [p.figdir,'figure2AB.pdf'];
pathtitle(fname)
fprintf('%s\n',savprint(h,fname));
close all

%%