function [h,scrsz]=ffigure
% ffigure - build fullscreen-size figure window
% function [h,scrsz]=ffigure
% 
% 
% [Outputs]
%     h: figure handle
%     scrsz: screen size
%     
%     
% Tomoyasu Horikawa horikawa-t@atr.jp
% 
% 
scrsz = get(0,'ScreenSize');
h=figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]);

