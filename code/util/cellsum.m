function sumofX=cellsum(X)
% cellsum -- sum all the cell component 
% sumofX=cellsum(X)
% 
% [Input]
%   -X : cell 
% 
% 
% 
% [Output]
%   -sumofX: sum of Input
% 
% 
% 
% 
% 
% 
% Created by Tomoyasu Horikawa 2009/12/16
% 

sumofX=zeros(size(X{1}));
for itr=1:numel(X)
    sumofX=sumofX+X{itr};
end
