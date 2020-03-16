function [slope,intercept,stats]=demingRegression(x,y,intercept0,lambda,null_slope,tail,nfold)
% [slope,intercept,stats]=demingRegression(x,y,intercept0,lambda,null_slope,tail,nfold)
% Perform deming regression
%
% [input]
%  x: input variable 1
%  y: input variable 2
%  intercept0: if 1 the intercept is to be constrained to zero (default = 1)
%  lambda: ratio of the variance of y and x (default = 1)
%  null_slope: slope for nulll hypothesis (default = 1)
%  tail: tail for significance test 'both' or 'left' or 'right' (default = 'both')
%  nfold: the number of fold for  (default = # of samples)
%
%
% [Output]
%  slope: slope
%  intercept: intercept
%  stats: statistics for significance test
%   .slopes: sloples estimated by the jack-knife method
%   .se: standard error of the slope estimated by the jack-knife method
%   .tval: t-values estimated by the jack-knife method
%   .pval: p-values estimated by the jack-knife method
%
% Tomoyasu Horikawa horikawa-t@atr.jp 20181219
%

%% setting
if ~exist('intercept0','var') || isempty(intercept0)
    intercept0 = 1;
end

if ~exist('lambda','var') || isempty(lambda)
    lambda = 1;
end

if ~exist('nfold','var') || isempty(nfold)
    nfold = length(x);
end
if nfold < 2
    error('nfold must be larger than 2.')
end
if nfold > length(x)
    warning('nfold must be shorter than x.')
    nfold = length(x);
end

%% estimate
if intercept0 == 1
    sxx = sum(x.^2);
    syy = sum(y.^2);
    sxy = sum(x.*y);
    slope = (-lambda*sxx+syy+sqrt(4*lambda*sxy^2+(lambda*sxx-syy)^2))/(2*sxy);
    intercept = 0;
else
    c_xy = cov(x,y);
    sxx = c_xy(1,1);
    sxy = c_xy(1,2);
    syy = c_xy(2,2);
    slope = (-lambda*sxx+syy+sqrt(4*lambda*sxy^2+(lambda*sxx-syy)^2))/(2*sxy);
    intercept = mean(y)-slope*mean(x);
end

% calculate se and sloples by the jack-knife method for significance test
% of slope over null_slope.
% It can take a few minutes, if sample size is large.
if nargout > 2
    if ~exist('null_slope','var') || isempty(null_slope)
        null_slope = 1;
    end
    if ~exist('tail','var') || isempty(tail)
        tail = 'both';
    end
    
    % jackknife estimates 
    slopes = zeros(nfold,1);
    n = length(x);
    if n == nfold
        held_outs = [1:n;1:n];
    else
        steps = round(linspace(1,n,nfold+1));
        held_outs = [steps(1:end-1);[steps(2:end-1)-1,steps(end)]];
    end
    allIdx = 1:n;
    for jk = 1:nfold
        idx = setdiff(allIdx,held_outs(1,jk):held_outs(2,jk));
        slopes(jk) = demingRegression(x(idx),y(idx),intercept0,lambda);
    end
    
    se = sqrt((nfold-1)*mean((slopes-slope).^2));
    tval = (slope-null_slope)/se;
    switch tail
        case 'both' % 'both or 'ne'
            pval = 2*tcdf(-abs(tval),nfold-1);
        case 'right' % 'right' or 'gt'
            pval = tcdf(-tval,nfold-1);
        case 'left' % 'left' or 'lt'
            pval = tcdf(tval,nfold-1);
    end
    
    % summary
    stats.slopes = slopes;
    stats.se = se;
    stats.tval = tval;
    stats.pval = pval;
end

%%

