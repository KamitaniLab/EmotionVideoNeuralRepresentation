function [nrow,ncol,ord]=setrc2(nDraw,direction,rc,off,centeroff)
% setrc- set the nrow and ncol for square
%
%
% [Inputs]
%     -nDraw: number of drawing
%     -direction: direction of order default='ltr'(left-top-right)
%      'ltr': left-top-right
%      'ltd': left-top-down
%      'rtl': right-top-left
%      'rtd': right-top-down
%      'lbr': left-bottom-right
%      'lbu': left-bottom-up
%      'rbl': right-bottom-left
%      'rbu': right-bottom-up
% [Optional inputs]
%     if r and c is pre-specified, only the order was calculated
%     -rc: specified row and column
%     -off: remove outframe: n of remove
%           [r,c]
%     -centeroff: remove center: n of remove
%          remove, 1:column, 2:row, 3:both, 4:diagonal, 5:upper triangle
%
% [Outputs]
%     -nrow:number of rows
%     -ncol:number of colmns
%     -ord: order for drawing subplot
%
%
%     [r,c,o]=setrc(12,'ltr')
%     [r,c,o]=setrc(12,'ltd')
%     [r,c,o]=setrc(12,'rtl')
%     [r,c,o]=setrc(12,'rtd')
%     [r,c,o]=setrc(12,'lbr')
%     [r,c,o]=setrc(12,'lbu')
%     [r,c,o]=setrc(12,'rbl')
%     [r,c,o]=setrc(12,'rbu')
%
%
% Tomoyasu Horikawa horikawa-t@atr.jp 2010/7/28
% modified by Tomoyasu Horikawa horikawa-t@atr.jp 2015/12/18
%

if ~exist('direction','var') || isempty(direction)
    direction='ltr';
end
if ~exist('rc','var') || isempty(rc)
    nrow=floor(nDraw/ceil(sqrt(nDraw)));
    ncol=ceil(nDraw/floor(nDraw/ceil(sqrt(nDraw))));
else
    nrow=rc(1);
    ncol=rc(2);
end
mat=makeM(1:(nrow*ncol),[nrow,ncol],2);
if exist('off','var') && ~isempty(off) && any(off~=0)
    if length(off)==1
        for i=1:off(1)
            mat(1,:)=[];
            mat(end,:)=[];
            mat(:,1)=[];
            mat(:,end)=[];
        end
    else
        for i=1:off(1)
            mat(1,:)=[];
            mat(end,:)=[];
        end
        for i=1:off(2)
            mat(:,1)=[];
            mat(:,end)=[];
        end
    end
end
if exist('centeroff','var') && ~isempty(centeroff) && centeroff~=0
    [cr,cc]=size(mat);
    switch centeroff
        case 1
            mat(:,ceil(cc/2))=[];
        case 2
            mat(ceil(cr/2),:)=[];
        case 3
            mat(ceil(cr/2),:)=[];
            mat(:,ceil(cc/2))=[];
        case 4
            for i=1:size(mat,1)
                mat(i,i)=NaN;
            end
        case 5
            for i=1:size(mat,1)
            for j=1:size(mat,1)
                if i<=j
                mat(i,j)=NaN;
                end
            end
            end
    end
end
switch direction
    case 'ltr'% left-top-right
        ord=asvector(mat,2);
    case 'ltd'% left-top-down
        ord=asvector(mat,1)';
    case 'rtl'% right-top-left
        matx=mat(:,end:-1:1);
        ord=asvector(matx,2);
    case 'rtd'% right-top-down
        matx=mat(:,end:-1:1);
        ord=asvector(matx,1)';
    case 'lbr'% left-bottom-right
        matx=mat(end:-1:1,:);
        ord=asvector(matx,2);
    case 'lbu'% left-bottom-up
        matx=mat(end:-1:1,:);
        ord=asvector(matx,1)';
    case 'rbl'% right-bottom-left
        matx=mat(end:-1:1,end:-1:1);
        ord=asvector(matx,2);
    case 'rbu'% right-bottom-up
        matx=mat(end:-1:1,end:-1:1);
        ord=asvector(matx,1)';
end
ord=ord(~isnan(ord));

% nrow=floor(nDraw/round(sqrt(nDraw)));
% ncol=ceil(nDraw/floor(nDraw/round(sqrt(nDraw))));

% nrow=round(nDraw/round(sqrt(nDraw)));
% ncol=round(nDraw/round(nDraw/round(sqrt(nDraw))));

% fprintf('subplot(nrow,ncol,itr)\n')