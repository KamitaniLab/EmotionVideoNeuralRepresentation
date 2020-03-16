function V=vecInc(input, incre)
% vecInc--increase the # of each components of input by each components of incre 
% 
%  e.g,
%   V=vecInc(1:3,1:3)
%   V =
%      1     2     2     3     3     3
% 
% 
% related function:
%   merge,
% created by HORIKAWA tomoyasu 09/09/03
V=cell(length(input),1);
if numel(incre)==1
    incre=repmat(incre,length(input),1);
end
for i=1:length(input)
    V{i}=repmat(input(i),1,incre(i));
end
V=merge(V,2);

