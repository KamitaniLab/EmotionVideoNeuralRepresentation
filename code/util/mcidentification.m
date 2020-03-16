function [cr,ci,pw_wins] = mcidentification(simmat,labels,ncandidates,nrepetitions)
% mcidentification -- perform multiclass identification from similarity matrix
% function mcidentification(simmat,labels,ncandidates)
%
% [Inputs]
%     -simmat: Similarity matrix [pred x correct]
%     -labels: index matrix
%     -ncandidates: numbers of candidates [2~ncorrect]
%     -nrepetitions: numbers of repetitions for multi-class identification
%
% [Outputs]
%     -cr: all correct rate
%     -ci: confidence inter vals
%     -pw_wins: winnter index for pairwise analysis
%
%
% Written by Tomoyasu Horikawa horikawa-t@atr.jp 2019/3/5
%

%% settings

if ~exist('ncandidates','var') || isempty(ncandidates)
    fprintf('No candidnate info. Do pairwise identificaion.\n')
    ncandidates = 2;
end
if ~exist('nrepetitions','var') || isempty(nrepetitions)
    nrepetitions = 200;
end

% remove nan if any
if any(any(isnan(simmat),1))
    nanIdx = isnan(simmat);
    simmat(any(nanIdx,2),:) = [];
    labels(any(nanIdx,2)) = [];
end

nSamp = size(simmat,1);
candidateIdx = 1:size(simmat,2);
% calculate order
[sorted,ord] = sort(simmat,2,'descend');

% initialize
cr = cell(length(ncandidates),1);
ci = cell(length(ncandidates),1);
for ix = 1:length(ncandidates)
    cr{ix} = zeros(nSamp,1);
    ci{ix} = zeros(nSamp,1);
end
pw_wins = zeros(nSamp,1);

%% identification
for i = 1:nSamp
    similars = ord(i,1:find(ord(i,:) == labels(i))); % find the more similar candidates
    candidate = candidateIdx(candidateIdx ~= labels(i)); % exclude true candidate
    
    pw_wins(i,candidateIdx ~= labels(i)) =~ ismember(candidate,similars(1:end-1));
    evalIdx = find(candidateIdx ~= labels(i));
    for ix = 1:length(ncandidates)
        if ncandidates(ix) == 2
            [ci{ix}(i,:),cr{ix}(i,:)] = ciestim3(pw_wins(i,evalIdx),2); % correct rate
        else
            cr_tmp = zeros(nrepetitions,1);
            for ixx = 1:nrepetitions
                rand('seed',ix*ixx)
                randIdx = randsample(evalIdx,ncandidates(ix)-1);
                cr_tmp(ixx) = all(pw_wins(i,randIdx) == 1); % correct rate
            end
            [ci{ix}(i),cr{ix}(i)] = ciestim3(cr_tmp,1); % correct rate
        end
    end
end


%%

