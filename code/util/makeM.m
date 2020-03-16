function B=makeM(A,rowcol,direction)
% makeM -- builds the matrix which consist of component of vecor A
% function B=makeM(A,rowcol,direction)
% Input
%   A: vector.
%   rowcol: two components vector[row col]
%           size(B) where B is matrix can be directly applied
%   direction: direction of setting component. direction ==1 «@direction ==2 -> (default=1)
% The components are entered from left to right, not top to bottom
% e.g.,
%  B=makeM(1:10,[2,5],2)
%   B =
%      1     2     3     4     5
%      6     7     8     9    10
%   A=makeM(1:10,size(B),2)
%   A =
%      1     2     3     4     5
%      6     7     8     9    10
%   A=makeM(1:10,size(B),1)
%      1     3     5     7     9
%      2     4     6     8    10
%
%
% created by  HORIKAWA tomoyasu 09/09/03
% modified by HORIKAWA tomoyasu 09/09/03
% modified by HORIKAWA tomoyasu 09/10/23 switch the direction 1<->2
switch direction
    case 1
        for i=1:rowcol(2)
            B(:,i)=A((1:rowcol(1))+(i-1)*rowcol(1));
        end
    case 2
        for i=1:rowcol(1)
            B(i,:)=A((1:rowcol(2))+(i-1)*rowcol(2));
        end
    otherwise
end



