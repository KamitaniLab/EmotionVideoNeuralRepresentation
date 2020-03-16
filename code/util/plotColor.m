function c=plotColor(ind,n,type,inv)
% plotColor -- create color matrix of type colors
% function c=plotColor(n,type)
%
% [Outputs]
%   c:color color
%
%
%
% Tomoyasu Horikawa horikawa-t@atr.jp 2011/12/4
%
if ~exist('n','var') || isempty(n)
    n=10;
end
if ~exist('type','var') || isempty(type)
    type='hot';
end
if ~exist('inv','var') || isempty(inv)
   inv=0;
end

ignoreRange=0; % to avoid white color
mxind=64;
rang=mxind-ignoreRange;
onset=1;
switch type
    case {'Paired_r','rainbow_r','Spectral_r'}
        cmat = colmap(type(1:end-2));
        cmat = cmat(end:-1:1,:,:);
    case {'Paired','viridis','terrain','rainbow','Spectral','Pastel1_r','RdBu_r','coolwarm','RdGy_r'}
        cmat = colmap(type);
    case 'spring'
        cmat=spring;
    case 'summer'
        cmat=summer;
    case 'autumn'
        cmat=autumn;
    case 'winter'
        cmat=winter;
    case 'cool'
        cmat=cool;
    case 'hot'
        cmat=hot;
        ignoreRange=10;
        rang=mxind-ignoreRange;
    case 'hsv'
        cmat=hsv;
        cmat=cmat(end:-1:1,:);
        onset=10;
    case 'bone'
        cmat=bone;
        ignoreRange=10;
        rang=mxind-ignoreRange;
    case 'pink'
        cmat=pink;
        ignoreRange=10;
        rang=mxind-ignoreRange;
    case 'copper'
        cmat=copper;
    case 'copbon'
        tmp=bone;
        cmat=[tmp(end-10:-1:10,:);copper];
        mxind=size(cmat,1);
        rang=mxind-ignoreRange;
    case 'rbb'
        cmat=rbb;
    case 'rwb'
        cmat=rwb;
    case 'ybm'
        cmat=ybm;
    case {'jet','jetpa'}
        cmat=jet;
    otherwise
        cmat=hsv;
end
c=cell(n,1);
sampling=round(linspace(onset,rang,n));
for itr=1:n
    c{itr,1}=cmat(sampling(itr),:);
end
switch type
    case {'jetpa'}
        seq = 1:3;
        for itr=1:n
            [mx,mxord] = max(c{itr,1});
            c{itr,1}(seq(seq~=mxord)) = c{itr,1}(seq(seq~=mxord));
            c{itr,1}(seq(seq==mxord)) = c{itr,1}(seq(seq==mxord))*0.5;
        end
end

switch type
    case 'rbow'
        c={...
            [1,0,0]     % red
            [1,.6,0]    % orange
            [1,1,0]     % yerllow
            [.8,1,.4]   % light green
            [.2,1,.6]   % emerald
            [0,1,0]     % green
            [0,.6,0]    % dark green
            [.3,.5,1]   % light blue
            [0,1,1]     % cyan
            [0,0,.9]    % blue
            [.6,.2,.8]  % purple
            [1,0,1]     % magenda
            [1,.6,1]    % pink
            [0,0,0]     % black
            [.3,.3,.3]  % gray1
            [.6,.6,.6]  % gray2
            [.9,.9,.9]  % gray3
            [1,0,0]     % red
            [1,.6,0]    % orange
            [1,1,0]     % yerllow
            [.8,1,.4]   % light green
            [.2,1,.6]   % emerald
            [0,1,0]     % green
            [0,.6,0]    % dark green
            [.3,.5,1]   % light blue
            [0,1,1]     % cyan
            [0,0,.9]    % blue
            [.6,.2,.8]  % purple
            [1,0,1]     % magenda
            [1,.6,1]    % pink
            [0,0,0]     % black
            [.3,.3,.3]  % gray1
            [.6,.6,.6]  % gray2
            [.9,.9,.9]  % gray3
            };
        mxind=17;
        c=c(end:-1:1);
    case 'rbow10'
        c={...
            [1,0,0]     % red
            [1,.6,0]    % orange
            [.8,1,.4]   % light green
            [0,.6,0]    % dark green
            [0,1,1]     % cyan
            [0,0,.9]    % blue
            [.6,.2,.8]  % purple
            [.1,.1,.1]  % gray1
            [.4,.4,.4]  % gray2
            [.7,.7,.7]  % gray3
            };
        mxind=10;
        c=c(end:-1:1);
    case {'rbsemo'} %  red blue and sky for emotion
        c{1}=[235, 1, 22]/255; %
        c{2}=[18, 94, 240]/255; %
        c{3}=[10, 168, 247]/255; %
    case {'rbemo'} %  red and blue for emotion
        c{1}=[235, 1, 22]/255; %
        c{2}=[18, 94, 240]/255; %
    case {'pb'} % sky pink blue
        c{2}=[250,120,210]/255; %
        c{1}=[39,169,225]/255; %
    case {'spb'} % sky pink blue
        c{1}=[197,234,245]/255; % sky
        c{3}=[250,120,210]/255; % pink
        c{2}=[39,169,225]/255; % blue
    case {'pa','pastel'}
        c{1}=[1,0.7,0.7];
        c{2}=[0.7,1,1];
        c{3}=[0.7,0.9,0.7];
        c{4}=[0.4,0.4,0.4];
        c{5}=[0.7,0.7,1];
    case {'pa2','pastel2'}
        c{1}=[1,0.8,0.5];
        c{2}=[1,0.7,1];
        c{3}=[0.7,0.9,0.7];
        c{4}=[0.4,0.4,0.4];
        c{5}=[0.7,0.7,1];
    case {'pa3','pastel3'}
        c{1}=[1,0.7,0.7];
        c{2}=[0.7,1,1];
        c{3}=[0.7,0.9,0.7];
        c{4}=[0.4,0.4,0.4];
        c{5}=[0.7,0.7,1];
        c{6}=[1,0.8,0.5];
        c{7}=[1,0.7,1];
    case {'pa4','pastel4'}
        c{1}=[1,0.7,0.7];
        c{2}=[0.7,1,1];
        c{3}=[0.7,0.9,0.7];
        c{4}=[0.7,0.7,1];
    case {'pa5','pastel5'}
        c{1}=[1,0.7,0.7];
        c{2}=[0.7,1,1];
        c{3}=[0.7,0.9,0.7];
        c{4}=[0.7,0.7,1];
        c{5}=[0.4,0.4,0.4];
        c{6}=[1,0.8,0.5];
        c{7}=[1,0.7,1];
        c{8}=[0.6,0.6,0.6];
        c{9}=[0.8,0.8,0.8];
    case {'pa6','pastel6'}
        c{1}=[1,0.8,0.5];
        c{2}=[1,0.7,1];
        c{3}=[0.7,0.9,0.7];
        c{4}=[1,0.7,0.7];
        c{5}=[0.7,0.7,1];
    case {'pa0'}
        c{1}=[0.4,0.4,0.4]; % dark gray
        c{2}=[1,0.7,0.7]; % red
        c{3}=[0.7,1,0.7]; % green
        c{4}=[0.9,0.5,1]; % purple
    case {'suteki4','s4'}
        c{1}=[0.1,0.1,0.1]; % dark gray
        c{2}=[1,0.2,0.2]; % red
        c{3}=[0.2,1,0.2]; % green
        c{4}=[0.7,0.2,1]; % purple
    case {'suteki7','s7'}
        c{1}=[0.1,0.1,0.1]; % dark gray
        c{2}=[1,0.2,0.2]; % red
        c{3}=[0.2,1,0.2]; % green
        c{4}=[0.2,0.2,1]; % blue
        c{5}=[1,0.7,0.2]; % orange
        c{6}=[0.7,0.2,1]; % purple
        c{7}=[0.5,1,1]; % cyan
    case {'orb'}
        c{3}=[255,150,0]/255; %
        c{2}=[255,0,0]/255; %
        c{1}=[0,0,255]/255; %
end
c=c(end:-1:1);
if nargin>0 && ~isempty(ind)
%     c=c{mod(ind-1,mxind)+1};
    c=c{mod(ind-1,size(c,1))+1};% TH 20140306
end
if inv
    c=c(end:-1:1);
end


%% Example
%{
close all
h=figure;
% col=rbow;
nline=10;
nperiods=20;
col={'spring','summer','autumn','winter','cool','hot','hsv','bone','pink','copper','copbon','rbb','jet','rbow','ybm',...
    'Paired','viridis','terrain','rainbow','Spectral','Pastel1_r','RdBu_r','coolwarm'};
ncol=length(col);
[r,c]=setrc(length(col));
for citr=1:ncol
    subplot(r,c,citr)
    for itr=1:nline
        hold on
        plot(randn(1,nperiods)-itr*3,'Color',plotColor(itr,nline,col{citr}),'LineWidth',3)
    end
    axis off image
    title(col{citr})
end
ffine(h)
setdir('/home/mu/horikawa-t/util/god/figure/');
savprint(h,'/home/mu/horikawa-t/util/god/figure/plotColor.pdf')
%}



