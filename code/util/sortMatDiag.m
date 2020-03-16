function [sortedMat,rowOrd,colOrd] = sortMatDiag(d,discoutrate,rowSortType)
% function [sortedMat,rowOrd,colOrd] = sortMatDiag(d)
% % sort matrix such that digonal lines tended to have large values
%
% [input]
%   d: original matrix
%   discoutrate: discout rate for judging to move next row. default = 1;
%   rowSortType: how to sort rows['max', 'entropy']. default = 'max'
%
% [output]
%   sortedMat: sorted matrix
%   rowOrd: row order for sorting
%   colOrd: column order for sorting
%
% written by horikawa-t@atr.jp 20190827
%
%% parameters
if ~exist('discoutrate','var') || isempty(discoutrate)
    discoutrate = 1; % discount rate
end

if ~exist('rowSortType','var') || isempty(rowSortType)
    rowSortType = 'max'; % discount rate
end

%%
co0 = 1:size(d,2);
ro0 = 1:size(d,1);
corest = 1:size(d,2);
rorest = 1:size(d,1);
co = 1:size(d,2);
ro = 1:size(d,1);
cox = [];
rox = [];
ncolassigned = 0;
d_org = d;
switch rowSortType
    case 'colmax'
        for i = 1:size(d,1)
            indMat = d(i,:) == max(d(i,:));
            cidx = find(sum(indMat,1),1,'first');
            co_tmp = [corest(cidx),setxor(corest(cidx),corest)];
            
            drest = d_org((i+1):end,co_tmp(2:end));
            if isempty(drest)
                cox = [cox,corest];
                break
            end
            
            corest = setxor([cox,corest(i+1:end)],co);
            d = d_org((i+1):end,corest);
            if isempty(d)
                cox = [cox,corest];
                rox = 1:size(d_org,1);
                break
            end
        end
        
    otherwise
        for i = 1:size(d,1)
            indMat = d == max(d(:));
            cidx = find(sum(indMat,1),1,'first');
            switch rowSortType
                case 'max'
                    ridx = find(indMat(:,cidx),1,'first');
                case 'entropy'
                    h = zeros(size(d,1),1);
                    for ix = 1:size(d,1)
                        h(ix) = entropy0(d(ix,:)./sum(d(ix,:)));
                    end
                    [mi,ord] = sort(h,'ascend');
                    ridx = ord(1);
            end
            co_tmp = [corest(cidx),setxor(corest(cidx),corest)];
            ro_tmp = [rorest(ridx),setxor(rorest(ridx),rorest)];
            
            drest = d_org(ro_tmp(2:end),co_tmp(2:end));
            if isempty(drest)
                cox = [cox,corest];
                rox = [rox,rorest];
                break
            end
            [sorted,ord] = sort(d_org(ro_tmp(1),corest),'descend');
            useIdx = find(sorted >= max(drest(:))*discoutrate);
            
            ncolassigned = ncolassigned + length(useIdx);
            
            co = [cox,corest(ord(useIdx)),setxor([cox,corest(ord(useIdx))],co)];
            ro = [rox,rorest(ridx),setxor([rox,rorest(ridx)],ro)];
            
            cox = [cox,corest(ord(useIdx))];
            rox = [rox,rorest(ridx)];
            
            corest = setxor([cox,corest(ord(useIdx))],co);
            rorest = setxor([rox,rorest(ridx)],ro);
            d = d_org(rorest,corest);
            
            if isempty(d)
                cox = [cox,corest];
                rox = [rox,rorest];
                break
            end
        end
end
sortedMat = d_org(rox,cox);
rowOrd = rox;
colOrd = cox;

