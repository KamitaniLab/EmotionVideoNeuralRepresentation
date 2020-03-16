function h = subplottight(rows, cols, ax, margin)
% subplottight    Subplot but it removes the spacing between axes
%
% This file is a part of BrainDecoderToolbox2.
%
% Inputs:
%
% - rows   : number of rows
% - cols   : number of columns
% - ax     : requested subplot axis
% - margin : amount of margin (in percent) separating axes
%            scalar or [xmargin, ymargin]
%
% Outputs:
%
% - h : handle to the axis
%
% Note:
%
% This function is based on `subplottight` (author unknown).
%

if nargin < 4
    xmargin = 0.01;
    ymargin = 0.01;
else
    if length(margin) == 2
        xmargin = margin(1);
        ymargin = margin(2);
    else
        xmargin = margin;
        ymargin = margin;
    end
end

ax = ax-1;
x = mod(ax,cols)/cols;
y = (rows-fix(ax/cols)-1)/rows;

h = axes('Units', 'normalized', ...
         'Position', [x + xmargin/cols, ...
                      y + ymargin/rows, ...
                      1/cols - 2*xmargin/cols, ...
                      1/rows - 2*ymargin/rows]);
