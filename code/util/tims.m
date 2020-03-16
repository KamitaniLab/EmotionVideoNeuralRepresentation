function s=tims(s)
%% tims -- outputs elapsed time from las tic
%
% [Inputs]
%     -s: if 1 sprintf is used. if 0 fprintf is used. (default=0)
% [Outputs]
%     -s: 0 or sprint text
%
% Tomoyasu Horikawa 20130730
%
try
    if exist('s','var') && s==1
        s=sprintf('[%2dd %2dh %2dm %2ds]:',floor(toc/86400),mod(floor(toc/3600),24),mod(floor(toc/60),60),floor(mod(toc,60)));
    else
        fprintf('[%2dd %2dh %2dm %2ds]:\n',floor(toc/86400),mod(floor(toc/3600),24),mod(floor(toc/60),60),floor(mod(toc,60)));
    end
catch me
    fprintf('tic is not set.\n')
    tic
    try
        fprintf('NOW:%s\n',datestr(datevec(now)))
    end
end