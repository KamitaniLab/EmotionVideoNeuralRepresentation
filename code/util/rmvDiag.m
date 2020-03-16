function mat=rmvDiag(sqMat,direction)
% rmvDiag - remove diagonal component from square matrix
% function rmvDiag(sqMat,direction)
%
% [Input]
%  sqMat: square matrix
%  direction:; direction to move the rest (default = 'x' or 1)
%
%
% [Output]
%  mat: matrix
%
% witten by Tomoyas Horikawa horikawa.t@gmail.com 2016/03/04

if ~exist('direction','var') || isempty(direction)
    direction = 'x';
end
switch direction
    case {1,'x'}
        if iscell(sqMat)
            mat= cell(size(sqMat,1),size(sqMat,2)-1);
        else
            mat= zeros(size(sqMat,1),size(sqMat,2)-1);
        end
        seq=1:size(sqMat,2);
        for i=1:size(sqMat,1)
            mat(i,:)=sqMat(i,seq~=i);
        end
        
    case {2,'y'}
        if iscell(sqMat)
            mat= cell(size(sqMat,1)-1,size(sqMat,2));
        else
            mat= zeros(size(sqMat,1)-1,size(sqMat,2));
        end
        seq=1:size(sqMat,1);
        for i=1:size(sqMat,2)
            mat(:,i)=sqMat(seq~=i,i);
        end
        
end








