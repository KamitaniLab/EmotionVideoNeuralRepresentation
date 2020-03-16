function ffine(h)
% ffine -- set paramerters of figure
% 
% [Inputs]
%     -h:handle
%     
%     
%     
%     
% [Ex.]
% close all
% h=figure;
% n=10;
% x=1:n;
% y1=randn(n,1);
% y2=randn(n,1);
% plot( x, y1 );
% hold on;
% plot( x, y2 );
% ffine(h)
%     
%     
%     
% Written By Tomoyasu Horikawa horikawa-t@atr.jp 2011/08/27
% 
%%
fh = figure(h); % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white
% set(fh, 'color', 'black'); % sets the color to white
set(gca, 'Box', 'off' ); % here gca means get current axis
set(gca, 'TickDir', 'out');
legend BOXOFF
% set(h1, 'LineStyle', '-', 'LineWidth', 1.0, 'Color', 'Black');
% set(h1, 'Marker', 'o', 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 8.0);
% set(h2, 'LineStyle', '?', 'LineWidth', 1.0, 'Color',' Black');
% set(h2, 'Marker', 's', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', [0 0 0], 'MarkerSize', 8.0);
% set(gca, 'XTick', [1:10], 'YTick', [0 10 20 30]);



%%