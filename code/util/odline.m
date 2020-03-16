function odline(ax)
% odline -- draw off diagonal line
%
% [Input]
%   -ax:axis
%
%
%
%
%
%
%
%
% Written by Tomoyasu horikawa horikawa-t@atr.jp 2011/10/05
%
%

hld=ishold;
hold on
min2max=[min([ax(1),ax(3)]),max([ax(2),ax(4)])]+...
    ([-abs(min([ax(1),ax(3)])),abs(max([ax(2),ax(4)]))])*0.1;
h=plot(min2max,min2max,'-k');
x=get(h,'Annotation');
axis([min2max,min2max])
if iscell(x)
    for i=1:length(x)
        x{i}.LegendInformation.IconDisplayStyle='off';
    end
else
    x.LegendInformation.IconDisplayStyle='off';
end

if ~hld
    hold off
end

