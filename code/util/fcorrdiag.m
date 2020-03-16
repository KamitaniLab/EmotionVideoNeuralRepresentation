function [cor t p]=fcorrdiag(X,Y,tri)
% fcorr -- perform fast correlation calculation
% function cor=fcorr(X,Y)
%
% [Inputs]
%     -X: matrix or vector
%     -Y: matrix or vector
%     -tri:whther output unique triangle(default=0)
%
% [Outputs]
%     -cor:correlation between X and Y
%
%
%     Created By Tomoyasu Horikawa horikawa-t@atr.jp 2010/6/6
%
%

if ~exist('tri','var')
    tri=0;
end
N=size(X,1);

Xdim=size(X,2);
Ydim=size(Y,2);

if Xdim~=Ydim
    error('Data dimension must be the same\n')
end

cor=zeros(Xdim,1);
for i=1:Xdim
if ~exist('Y','var')
    cor(i)=cov(X(:,i))./(std(X(:,i))'*std(X(:,i)));
else
    if tri
        cor(i)=getTri(cov([X(:,i),Y(:,i)])./(std([X(:,i),Y(:,i)])'*std([X(:,i),Y(:,i)])));
    else
        %         cor=cov([X,Y])./...
        %             (std(X)'*std(Y));
        cor(i)=(((X(:,i)-repmat(mean(X(:,i)),N,1))'*(Y(:,i)-repmat(mean(Y(:,i)),N,1)))/(N-1))./...
            (std(X(:,i))'*std(Y(:,i)));
    end

end

if nargout>2
    % function p = pvalPearson(tail, rho, n)
    %PVALPEARSON Tail probability for Pearson's linear correlation.
    t = sign(cor) .* Inf;
    k = (abs(cor) < 1);
    t(k) = cor(k).*sqrt((N-2)./(1-cor(k).^2));
    tail='b';
    switch tail
        case 'b' % 'both or 'ne'
            p = 2*tcdf(-abs(t),N-2);
        case 'r' % 'right' or 'gt'
            p = tcdf(-t,N-2);
        case 'l' % 'left' or 'lt'
            p = tcdf(t,N-2);
    end
end
end
