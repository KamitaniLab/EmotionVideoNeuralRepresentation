function z=r2z(r)
% r2z -- fisher's z transform: correlation r is converted to z score
% z=r2z(r);
% 
% [Inputs]
%     -r:correlation coefficient
%     
% 
% 
% 
% [Outputs]
%     -z:zscore
% 
% 
% 
% Written by Tomoyasu Horikawa horikawa-t@atr.jp 2012/4/23
% 

z=log((1+r)./(1-r))/2;


