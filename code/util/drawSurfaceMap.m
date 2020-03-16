function coloredIm = drawSurfaceMap(scores,colMap,col,rang)
% function I = drawSurfaceMap(scores,colMap,color)
% - this code is used to make a surface map colored by score.
%
% [input]
%  scores: values used for drawing map
%  colMap: map information
%    .idxSurfMap: map for index
%    .colSurfMap: template map for output
%    .roiPosIdx: order index to match input and idxSurfMap
%    .roiSets: roiSets ordered in the idx in idxSurfMap
%  col: color schema: {'Paired','viridis','terrain','rainbow','Spectral','Spectral_r','Pastel1_r','RdBu_r','RdBu','hot','coolwarm','RdGy_r','RdGy'};
%  rang: range of colorbar, [min,max], if not provided, draw normalized map.
%
% [output]
%  coloredIm: colored surface map
%
%
% Note: currently, only HCP360 is suppored.
%
%
% written by horikawa-t 2020/03/14
%

%% make colored map
if ~exist('col','var')
    col = 'RdGy'; % 'RdBu_r'
end

nOnSurf = 360; % for HCP360

template = colMap.colSurfMap;
[h,w,d] = size(template);

% convert score to integer index
if ~exist('rang','var')
    % normalized
    scoreIdx = round(normal(scores)*(length(scores)-1)+1); % ranged 1:n elements
    % scoreIdx = calcrank(scores);
else
    scores(scores<=rang(1)) = rang(1);
    scores(scores>=rang(2)) = rang(2);
    scores = scores-rang(1);
    normscore = scores/(rang(2)-rang(1));
    scoreIdx = round(normscore*(length(scores)-1)+1); % ranged 1:n elements
end

% extend to n results
cm = colmap256(col);
cm = cm(round(linspace(1,size(cm,1),length(scoreIdx))),:);

cnt = 0;
% prepare lower blocks for out of surface results
nOver = length(scoreIdx)-nOnSurf;
cmat = cell(1,nOver);
if nOver > 0
    overWidth = round(w/nOver);
    short = w-overWidth*nOver;
    overWidths = repmat(overWidth,1,nOver);
    overWidths(end) = overWidths(end) + short;
    overHeight = round(h/10);
end
for i = 1:length(scoreIdx)
    if i > nOnSurf
        cnt = cnt+1;
        cmat{1,cnt}(:,:,1) = ones(overHeight,overWidths(cnt),1)*cm(scoreIdx(i),1)*255;
        cmat{1,cnt}(:,:,2) = ones(overHeight,overWidths(cnt),1)*cm(scoreIdx(i),2)*255;
        cmat{1,cnt}(:,:,3) = ones(overHeight,overWidths(cnt),1)*cm(scoreIdx(i),3)*255;
    else
        I1 = template(:,:,1);
        I2 = template(:,:,2);
        I3 = template(:,:,3);
        
        idx = colMap.idxSurfMap(:) == colMap.roiPosIdx(i);
        
        I1(idx) = cm(scoreIdx(i),1)*255;
        I2(idx) = cm(scoreIdx(i),2)*255;
        I3(idx) = cm(scoreIdx(i),3)*255;
        template(:,:,1) = I1;
        template(:,:,2) = I2;
        template(:,:,3) = I3;
    end
    %imagesc(template)
    %axis off image
    %drawnow
end

if nOver > 0
    coloredIm = [template;[cmat{:}]];
else
    coloredIm = template;
end

end
%%


