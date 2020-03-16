function draw_figureS3AB(p,encRes)
%
% This code is for drawing figure S3AB
%
%
%% settings
nSbj = length(p.sbjID);

%% draw figure S3A

%  draw results
close all
h = ffigure;
cnt = 0;

drawScores = {...
    'category'
    'categcontinuous'
    'dimension'
    'dim28binary'
    'dim28continuous'
    };


% display settings
r = 5;
c = 5;
[r,c,o] = setrc2(r*c,'ltr',[r,c],[1,1]);
fsize = 8;


acc = cell(length(drawScores),3);
axnam = cell(length(drawScores),1);
% prepare  results
for sbjitr = 1:length(p.sbjID)
    for ix = 1:length(drawScores)
        acc{ix,1} = encRes.(drawScores{ix}).pred_acc{sbjitr};
        axnam{ix} = sprintf('%s',drawScores{ix});
    end
    
    
    cnt = cnt+1;
    subplottight(r,c,o(cnt),0.15);
    violin(acc(:,1)','edgecolor','none','facecolor',[1,1,1]*0.4,'ci',0);
    set(gca,'FontSize',fsize)
    ylabel('Prediction accuracy (r)')
    ylim([-0.1,0.5])
    hline(0,'--k');
    axname(axnam)
    xticklabel_rotate([],45)
    title(sprintf('%s',p.sbjID{sbjitr}),'FontSize',fsize)
    ffine(h)
    
end
suptitle(sprintf('Comparison of various control features'));
fname = [p.figdir,'figureS3A.pdf'];
pathtitle(fname)
fprintf('%s\n',savprint(h,fname));

%% draw figure S3B
drawScores = {...
    {'category','categcontinuous'}
    {'category','dimension'}
    {'category','dim28binary'}
    {'category','dim28continuous'}
    };

% display settings
r = 5;
c = 5;
[r,c,o] = setrc2(r*c,'ltr',[r,c],[2]);
fsize = 8;
close all
h = ffigure;
cnt = 0;

% prepare  results
slopeDR = zeros(nSbj,length(drawScores));
angle = zeros(nSbj,length(drawScores));
scnames = cell(length(drawScores),1);
for ix = 1:length(drawScores)
    for sbjitr = 1:length(p.sbjID)
        ac1 = encRes.(drawScores{ix}{1}).pred_acc{sbjitr}';
        ac2 = encRes.(drawScores{ix}{2}).pred_acc{sbjitr}';
        [slopeDR(sbjitr,ix),int,st] = demingRegression(ac1,ac2,1,1,1,[],100);
        angle(sbjitr,ix) = rad2deg(atan(slopeDR(sbjitr,ix)));
        scnames{ix} = sprintf('%s',drawScores{ix}{2});
    end
end

cnt = cnt+1;
subplottight(r,c,o(cnt),0.05);
base = -0.03*3;
dif = 0.03;
for sbjitr = 1:length(p.sbjID)
    plot((1:length(drawScores))+base+dif*sbjitr,angle(sbjitr,:)-45,plotMarker(sbjitr),'Color',plotColor(sbjitr,nSbj,'pa'),'LineWidth',1);
    hold on
end
set(gca,'FontSize',fsize)
hline(0,'-k')
ylabel(sprintf('Deviation from equivalent angle'))
axname(scnames,1)
axname([-20:10:20],2,[-20:10:20])
xticklabel_rotate([],45)
text(1,15,'Ctrl','FontSize',fsize)
text(1,-15,'Category','FontSize',fsize)
legend(p.sbjID)
ylim([-20,20])
ffine(h)

suptitle(sprintf('Comparison of various control features'));
fname = [p.figdir,'figureS3B.pdf'];
pathtitle(fname)
fprintf('%s\n',savprint(h,fname));

%%
close all