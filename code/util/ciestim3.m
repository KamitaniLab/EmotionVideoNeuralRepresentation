function [hw,mu,sd,upper,lower,c]=ciestim3(x,rowcol,p,side,nan)
% ciestim -- estimate confidence interval
% function [hw,mu,sd,,upper,lower,c]=ciestim3(x,rowcol,p,side)
% 
% [Input]
%     -x:x vector[Nx1 or 1xN] or matrix[NxD]
%     -rowcol:direction of calculation default=1(row)
%     -p:confidence level (default = 0.95)
%     -side:onesided or twosided (default = onesided)
% 
% [Outputs]
%     -hw:half width of confidence interval
%     -mu:mean
%     -std:standard deviation
%     -upper:upper limit
%     -lower:lower limit
%     -c:cumulative probability
%     
% [Usage]
% 
% 
% [hw,u,l,c]=ciestim3(x,rowcol,p,side)
% 
% 
% [note]
% ref:http://kogolab.jp/elearn/hamburger/chap2/sec3.html
% 
% 
% Written by Tomoyasu Horikawa horikawa-t@atr.jp 2011/5/4
% 
% 

if ~exist('side','var') || isempty(side)
    side='onesided';
end
if ~exist('rowcol','var') || isempty(rowcol)
    rowcol=1;
end

if ~exist('p','var') || isempty(p)
    p=0.95;
end

if ~exist('nan','var') || isempty(nan)
    nan=1;
end

if nan
if rowcol==1
    nitr=size(x,2);
    n=zeros(1,nitr);
    for i=1:nitr
        n(i)=length(~isnan(x(:,i)));
    end
    c=zeros(1,nitr);
    hw=zeros(1,nitr);
else
    nitr=size(x,1);
    n=zeros(nitr,1);
    for i=1:nitr
        n(i)=length(~isnan(x(i,:)));
    end
    c=zeros(nitr,1);
    hw=zeros(nitr,1);
end
mu=nanmean(x,rowcol);
sd=nanstd(x,[],rowcol);
if 0 % TH170429
    for i=1:nitr
        switch side
            case 'onesided'
                c(i)=tinv(p,n(i)-1);
            case 'twosided'
                c(i)=tinv((1+p)/2,n(i)-1);
            otherwise
                error('Invalid ''side'' option.')
        end
        hw(i)=c(i).*sd(i)./sqrt(n(i));
    end
else
    switch side
        case 'onesided'
            c=tinv(p,n-1);
        case 'twosided'
            c=tinv((1+p)/2,n-1);
        otherwise
            error('Invalid ''side'' option.')
    end
    hw=c.*sd./sqrt(n);
end
upper=mu+hw;
lower=mu-hw;
else
n=size(x,rowcol);
switch side
    case 'onesided'
        c=tinv(p,n-1);
    case 'twosided'
        c=tinv((1+p)/2,n-1);
    otherwise
        error('Invalid ''side'' option.')
end
mu=mean(x,rowcol);
sd=std(x,[],rowcol);

hw=c.*sd./sqrt(n);
upper=mu+hw;
lower=mu-hw;
end

%%


