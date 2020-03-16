function m=plotMarker(ind,LM)
% plotMarker -- provide different marker
% function m=plotMarker(ind)
%
% [Inputs]
%   ind: index for differnt marker
%   LM: index for line/marker change order (default=1: line fixed marker change) 
%
% [Outputs]
%   m:marker
%
%
%
% Tomoyasu Horikawa horikawa-t@atr.jp 2012/11/28
%
if ~exist('LM','var') || isempty(LM)
    LM=1;
end
lineType={...
    '-'
    '-.'
    '--'
    ':'
    };
MarkerType={...
    'o'
    'd'
    '^'
    'v'
    's'
    'p'
    'h'
    '<'
    '>'
    '+'
    'x'
    '*'
    '.'
    ''
    };
nType=length(MarkerType)*length(lineType);
idx=mod(ind-1,nType)+1;

switch LM
    case 1
        line_idx=ceil(idx/length(MarkerType));
        marker_idx=mod(idx-1-(line_idx*length(MarkerType)),length(MarkerType))+1;
    case 2
        marker_idx=ceil(idx/length(lineType));
        line_idx=mod(idx-1-(marker_idx*length(lineType)),length(lineType))+1;
end
m=[MarkerType{marker_idx},lineType{line_idx}];

%% Example
%{
close all
h=ffigure;
nLine=56;
col='hsv';
for i=1:nLine
% plot(0:10,(0:10)*(i*10),plotMarker(i,1),'MarkerSize',10)
% plot(0:10,(0:10)*(i*10),plotMarker(i,2),'MarkerSize',10)
plot(0:10,(0:10)*(i*10),plotMarker(i,1),'Color',plotColor(i,nLine,col),'MarkerSize',10)
% plot(0:10,(0:10)*(i*10),plotMarker(i,2),'Color',plotColor(i,nLine,col),'MarkerSize',10)
drawnow
hold on
end
setdir('/home/mu/horikawa-t/util/god/figure/');
savprint(h,'/home/mu/horikawa-t/util/god/figure/plotMarker.pdf')
%}
%% Example 2
%{
close all
h=ffigure;
% col=rbow;
nperiods=10;
col={'spring','summer','autumn','winter','cool','hot','hsv','bone','pink','copper','copbon','rbb','jet','rbow'};
ncol=length(col);
[r,c]=setrc(length(col));
nline=56;
for citr=1:ncol
    subplot(r,c,citr)
    for itr=1:nline
        hold on
        plot(randn(1,nperiods)-itr*3,plotMarker(itr),'Color',plotColor(itr,nline,col{citr}),'LineWidth',1)
    end
    title(col{citr})
end
setdir('/home/mu/horikawa-t/util/god/figure/');
savprint(h,'/home/mu/horikawa-t/util/god/figure/plotColor2.pdf')
%}
%% Example 3
%{
close all
h=ffigure;
% col=rbow;
nperiods=10;
col={'hsv'};
ncol=length(col);
[r,c]=setrc(length(col));
nline=56;
for citr=1:ncol
    subplot(r,c,citr)
    for itr=1:nline
        hold on
        plot(randn(1,nperiods)-itr*3,plotMarker(itr),'Color',plotColor(itr,nline,col{citr}),'MarkerSize',10,'LineWidth',2)
    end
    title(col{citr})
    axis off
    ffine(h)
end
setdir('/home/mu/horikawa-t/util/god/figure/');
savprint(h,'/home/mu/horikawa-t/util/god/figure/plotColor3.pdf')
%}

%}


