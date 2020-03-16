function draw_figure4BC_S2B(p,encRes)
%
% This code is for drawing figure 4BC and S2B
%
%
%% settings

drawPairs = {...
    {'category','dimension'}
    {'category','vision'}
    {'category','semantic'}
    };
drawNames = {...
    {'cat','dim'}
    {'cat','vis'}
    {'cat','sem'}
    };
nSbj = length(p.sbjID);


%% draw figure 4B and S2B
% display settings
c = nSbj+2;
r = length(drawPairs) + 2;
[r,c,o] = setrc2(r*c,'ltd',[r,c],[1]);

% draw setting
close all
h = ffigure;
cnt = 0;
fsize = 8;
th = 0.111;
col = 'summer';

for sbjitr = 1:length(p.sbjID)
    % prepare  results
    acc = cell(1,length(p.scoreTypes));
    scnames = cell(length(p.scoreTypes),1);
    for ix = 1:length(p.scoreTypes)
        acc{ix} = encRes.(p.scoreTypes{ix}).pred_acc{sbjitr};
        scnames{ix} = p.scoreTypes{ix};
    end
    
    % draw feature decoding accuracy
    nVox = size(acc{ix},2);
    nmax = 2000; % 10000 voxels were shown in the paper
    
    for ixx = 1:length(drawPairs)
        cnt = cnt+1;
        subplottight(r,c,o(cnt),0.15);
        rand('seed',1)
        randidx = randsample(nVox,nmax);
        a1 = acc{ismember(scnames,drawPairs{ixx}{1})}(randidx)';
        a2 = acc{ismember(scnames,drawPairs{ixx}{2})}(randidx)';
        nnanidx = ~isnan(a1)&~isnan(a2);
        dscatter(a1(nnanidx),a2(nnanidx),'marker','.','msize',4);
        colormap(col)
        
        hold on
        [slopeDR,int,st] = demingRegression(acc{ismember(scnames,drawPairs{ixx}{1})}',acc{ismember(scnames,drawPairs{ixx}{2})}',1,1,1,[],100);
        plot(-0.2:0.01:0.6,(-0.2:0.01:0.6)*slopeDR,'-m')
        angle = rad2deg(atan(slopeDR));
        text(0.3,0.05,sprintf('a=%.1f',angle),'FontSize',fsize)
        %fprintf('%sVS%s:%s[p=%.5f]\n',drawNames{ixx}{1},drawNames{ixx}{2},p.sbjID{sbjitr},st.pval)
        text(0.3,-0.05,sprintf('p=%.5f',st.pval),'FontSize',fsize)
        set(gca,'FontSize',fsize)
        
        axis square
        xlabel(drawPairs{ixx}{1})
        ylabel(drawPairs{ixx}{2})
        xlim([-0.1,0.5])
        ylim([-0.1,0.5])
        odline(axis)
        hline(0,'-k');
        vline(0,'-k');
        hh = hline(th,'--k');
        set(hh,'Color',[1,1,1]*0.8)
        hh = vline(th,'--k');
        set(hh,'Color',[1,1,1]*0.8)
        title(sprintf('%s vs. %s: %s',drawNames{ixx}{1},drawNames{ixx}{2},p.sbjID{sbjitr}),'FontSize',fsize)
    end
    
end
suptitle(sprintf('Comparison of voxel encoding accuracy between feautre sets with slope'));
fname = [p.figdir,'figure4B_S2B.pdf'];
pathtitle(fname)
fprintf('%s\n',savprint(h,fname));

%% draw figure 4C

%  draw results
close all
h = ffigure;
cnt = 0;

% display settings
r = 5;
c = length(drawNames)+2;
[r,c,o] = setrc2(r*c,'ltr',[r,c],[2,1]);
fsize = 8;

acc = cell(length(p.sbjID),3);
axnam = cell(length(p.sbjID),1);

for ix = 1:length(drawPairs)
    for sbjitr = 1:length(p.sbjID)
        acc{sbjitr,1} = encRes.(drawPairs{ix}{1}).pred_acc{sbjitr};
        acc{sbjitr,2} = encRes.(drawPairs{ix}{2}).pred_acc{sbjitr};
        axnam{sbjitr} = sprintf('Subject %d',sbjitr);
    end
    
    cnt = cnt+1;
    subplottight(r,c,o(cnt),0.15);
    violin_h(acc(:,1)',acc(:,2)','edgecolor','none','facecolor',[1,1,1]*0.4,'ci',0);
    set(gca,'FontSize',fsize)
    ylabel('Prediction accuracy (r)')
    ylim([-0.1,0.5])
    hline(0,'--k');
    axname(axnam)
    xticklabel_rotate([],45)
    title(sprintf('%s vs.%s',drawNames{ix}{1},drawNames{ix}{2}),'FontSize',fsize)
    ffine(h)
    
end
suptitle(sprintf('Distribution of voxel encoding accuracy'));
fname = [p.figdir,'figure4C.pdf'];
pathtitle(fname)
fprintf('%s\n',savprint(h,fname));

%%
close all


