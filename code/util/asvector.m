function V=asvector(B,direction)
% asvector --  converts matrix to vector
% function V=asvector(B,direction)
% 
% [Inputs]
% B:matrix
% direction:direction ==1 ??direction ==2 -> (default=1)
%
% example:
% A={ 1 2 3
%     4 5 6 }
%   B=asvector(A,1)=[1 4 2 5 3 6]' modified by horikawa 090827 transpose
%   B=asvector(A,2)=[1 2 3 4 5 6]
% created by HORIKAWA tomoyasu 09/09/03
% modified by HORIKAWA tomoyasu 12/09/18
if exist('direction','var')==0
    direction=1;
end

switch direction
    case 1
        V=B(:);
    case 2
        b=B';
        V=b(:)';
    otherwise
        error('Invalid input.')
end

% modified 2012/9/18
% if iscell(B)
% if direction==1
%     V=cell(size(B,1)*size(B,2),1);
%     for n=1:(size(B,1))
%         for m=1:size(B,2)
%             V((m-1)*size(B,1)+n,1)=B(n,m);
%         end
%     end
% else
%     V=cell(1,size(B,1)*size(B,2));
%     for n=1:(size(B,2))
%         for m=1:size(B,1)
%             V(1,(m-1)*size(B,2)+n)=B(m,n);
%         end
%     end
% 
% end
% else    
% if direction==1
%     V=zeros(size(B,1)*size(B,2),1);
%     for n=1:(size(B,1))
%         for m=1:size(B,2)
%             V((m-1)*size(B,1)+n,1)=B(n,m);
%         end
%     end
% else
%     V=zeros(1,size(B,1)*size(B,2));
%     for n=1:(size(B,2))
%         for m=1:size(B,1)
%             V(1,(m-1)*size(B,2)+n)=B(m,n);
%         end
%     end
% 
% end
% end
