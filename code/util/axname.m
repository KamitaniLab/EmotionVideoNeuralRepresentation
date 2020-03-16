function axname(axnames,xy,raw,varargin)
% axname -- write x or y axis names
% function axname(axnames,xy,raw,varargin)
%
% [Inputs]
%     -axnames:axisnames
%     -xy:assign x axis==1 or y axis==2 (default==1,x axis)
%     -raw:whether the label index is straightforwardly used or not
%       --0 or range of axes
%
% Created By Tomoyasu Horikawa horikawa-t@atr.jp 2010/6/2
%
%
%

if ~exist('xy','var') || isempty(xy)
    xy=1;
end
if ~exist('raw','var') || isempty(raw)
    raw=0;
end

if raw==0
    switch xy
        case 2
            set(gca,'YTickLabel',axnames,'YTick',1:length(axnames),varargin{:})
            ylim([0.5,length(axnames)+0.5])
        otherwise
            set(gca,'XTickLabel',axnames,'XTick',1:length(axnames),varargin{:})
            xlim([0.5,length(axnames)+0.5])
    end
else
    switch xy
        case 2
            set(gca,'YTickLabel',axnames,'YTick',raw,varargin{:})
            %ylim([min(raw)-0.5,max(raw)+0.5])
        otherwise
            set(gca,'XTickLabel',axnames,'XTick',raw,varargin{:})
            %xlim([min(raw)-0.5,max(raw)+0.5])
    end
end
