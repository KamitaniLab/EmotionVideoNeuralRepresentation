function m=plotMarker2(ind)
% plotMarker2 -- provide different marker
% function m=plotMarker2(ind)
%
% [Inputs]
%   ind: index for differnt marker with small set no lines
%
% [Outputs]
%   m:marker
%
%
%
% Tomoyasu Horikawa horikawa-t@atr.jp 2020/2/24
%
lineType={...
    ''
    ''
    ''
    ''
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
    };
MarkerType={...
    'o'
    'd'
    '^'
    'v'
    's'
    '<'
    '>'
    '+'
    'x'
    'h'
    'p'
    '*'
    '<'
    '>'
    };
nType=length(MarkerType)*length(lineType);
idx=mod(ind-1,nType)+1;

line_idx=ceil(idx/length(MarkerType));
marker_idx=mod(idx-1-(line_idx*length(MarkerType)),length(MarkerType))+1;
m=[MarkerType{marker_idx},lineType{line_idx}];

%%

