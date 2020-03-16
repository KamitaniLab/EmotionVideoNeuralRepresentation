function [cor t p]=fcorr(X,Y,tri)
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
if ~exist('Y','var')
    cor=cov(X)./(std(X)'*std(X));
else
    if tri
        cor=getTri(cov([X,Y])./(std([X,Y])'*std([X,Y])));
    else
        %         cor=cov([X,Y])./...
        %             (std(X)'*std(Y));
        if size(X,2)==1 && size(Y,2)==1
            cor=(((X-mean(X))'*(Y-mean(Y)))/(N-1))./...
                (std(X)'*std(Y));
        else
            warning off
            mx = mean(X,1);
            my = mean(Y,1);
            cor=(((X-mx(ones(N,1),:))'*(Y-my(ones(N,1),:)))/(N-1))./...
                (std(X)'*std(Y));
            warning on
        end
    end
    
end

if nargout>2
    % function p = pvalPearson(tail, rho, n)
    % PVALPEARSON Tail probability for Pearson's linear correlation.
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
