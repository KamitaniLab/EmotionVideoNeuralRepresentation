function rank=calcrank(v,type)
% calcrank -- calculate rank
% 
% [Inputs]
%     v:vector
%     type:'ascend' or 'descend'=default
%     
% [Outputs]
%     rank:rank
% 
%     
%     
%     
%     
% Written by Tomoyasu Horikawa horikawa-t 2012/6/21    
%     
%%
if ~exist('type','var') || isempty(type)
    type='descend';
end
[sorted,ord]=sort(v,type);
[idx,rank]=ismember(1:length(v),ord);
