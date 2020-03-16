function p=pathset(p,onoff)
% pathset--set or remove the path 
%
% [Input]
%   p - path which you want set ot remove
%   onoff - 1: set the path , 0: remove the path
%
switch onoff
    case 1
        %genpath(p)
        addpath(genpath(p));
    case 0
        %genpath(p)
        rmpath(genpath(p));
    otherwise
        error('You should select 0 or 1 for onoff argument\n');
end
