function nX=normal(X,direction,criterion)
% normal - normalize the data
% nx=normal(X,direction,criterion)
%
% [Input]
%   -X: data
%   -direction: normalized dimension (default=1)
%    -1:mean of each column become 0
%    -2:mean of each row    become 0
%   -criterion: baseline (default='min')
%    set ... to 0    
%    -'min' : minimum of the data
%    -'mean': mean of the data
%    -'norm': make the vector norm 1
%
% [Output]
%   -nX: noramlized data
%
% [Related function]
%
% [Example]
%
% X=rand(10,10)*repmat(1:10,10,1)
% criterion='norm'
% subplot(3,1,1)
% plot(X)
% nXc=normal(X,1,criterion)
% nXr=normal(X,2,criterion)
% subplot(3,1,2)
% plot(nXc')
% subplot(3,1,3)
% plot(nXr)
%
% [note]
%
%
% Created  By: Tomoyasu Horikawa horikawa-t@atr.jp 2009/11/28

if ~exist('direction','var') || isempty(direction)
    if size(X,1)==1
        direction=2 ;
    else
        direction=1 ;
    end
end
if ~exist('criterion','var')
    criterion='mean' ;
    criterion='min' ;
end

switch criterion
    case 'mean'
        mu=mean(X,direction);
        s_d=std(X,[],direction);

        [nSamp,nDim]=size(X);
        [nrow,ncol]=size(mu);

        nX=(X-repmat(mu,nSamp/nrow,nDim/ncol))./...
            repmat(s_d,nSamp/nrow,nDim/ncol);
%        mean(nX,direction)
%        std(nX,[],direction)
    case 'min'
        ma=max(X,[],direction);
        mi=min(X,[],direction);

        [nSamp,nDim]=size(X);
        [nrow,ncol]=size(mi);

        nX=(X-repmat(mi,nSamp/nrow,nDim/ncol))./...
            repmat(ma-mi,nSamp/nrow,nDim/ncol);
%        max(nX,[],direction)
%        min(nX,[],direction)

    case 'norm'
        normX=sum(sqrt(X.*X),direction);
        
        [nSamp,nDim]=size(X);
        [nrow,ncol]=size(normX);

        nX=X./repmat(normX,nSamp/nrow,nDim/ncol);
%        sum(sqrt(nX.*nX),direction)

end

%%
