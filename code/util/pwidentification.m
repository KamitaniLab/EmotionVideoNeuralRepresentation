function [cr,ci,wins]=pwidentification(simmat,labels)
% pwidentification -- perform pairwise identification from similarity matrix
% matrix
% function pwidentification(simmat,labels,p)
% 
% [Inputs]
%     -simmat: Similarity matrix [pred x correct]
%     -labels: index matrix
% 
% 
% 
% [Outputs]
%     -cr: all correct rate
%     -ci: confidence inter vals
% 
% 
% Written by Tomoyasu Horikawa horikawa-t@atr.jp 2016/9/10
% 
% 
% 

%%

% remove nan if any
if any(any(isnan(simmat),1))
    nanIdx = isnan(simmat);
    simmat(any(nanIdx,2),:)=[];
    labels(any(nanIdx,2))=[];
end

nSamp=size(simmat,1);

candidateIdx=1:size(simmat,2);
% calculate order 
[sorted,ord]=sort(simmat,2,'descend');

cr=zeros(nSamp,1);
ci=zeros(nSamp,1);
wins=false(nSamp,1);
for i=1:nSamp
    similars=ord(i,1:find(ord(i,:)==labels(i))); % find the more similar candidates
    candidate=candidateIdx(candidateIdx~=labels(i)); % exclude true candidate
    wins(i,candidateIdx~=labels(i))=~ismember(candidate,similars(1:end-1));
    [ci(i,:),cr(i,:)]=ciestim3(wins(i,candidateIdx~=labels(i)),2); % correct rate
    
end


%%

%%

