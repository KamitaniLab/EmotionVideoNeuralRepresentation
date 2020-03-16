function draw_figure2DE(p,decRes)
%
% This code is for drawing figure 2DE
% In the paper, we used Pycortex (Gat et al., 2015) for drawing map.
% Here, similar visualization is mimicked without pycortex
%
%% settings

% prepare surface map
load([p.rootPath,'code/util/hcp360colMap.mat'],'hcp360colMap');

% ROI info.
hcpIdx = ~cellfun(@isempty,strfind(p.roiDescrip,'hcp180'));
subIdx = ismember(p.roiDescrip,{'thalamus','hyppocampus','hypothalamus','pallidum','brainstem','caudate','putamen','nuc_accum','amygdala','cerebellum'});
roiIdx = hcpIdx|subIdx;

nSbj = length(p.sbjID);


% map drawing settings
col = 'RdGy_r';%  col: color schema: {'Paired','viridis','terrain','rainbow','Spectral','Spectral_r','Pastel1_r','RdBu_r','RdBu','hot','coolwarm','RdGy_r','RdGy'};

scoreTypes = {'category','dimension'};

% get feature names
load([p.featdir,'category.mat']);
catname = L.featname';
load([p.featdir,'dimension.mat']);
dimname = L.featname';

%% draw figure 2DE
% figure settings
fsize = 8;

for scoritr = 1:length(scoreTypes)
    scoreType = scoreTypes{scoritr};
    switch scoreType
        case 'category'
            r = 9;
            c = 7;
            [r,c,o] = setrc2(r*c,'ltr',[r,c],[1]);
            res = decRes.category;
            featname = catname;
            figname = 'figure2D.pdf';
        case 'dimension'
            r = 9;
            c = 7;
            [r,c,o] = setrc2(r*c,'ltr',[r,c],[2]);
            res = decRes.dimension;
            featname = dimname;
            figname = 'figure2E.pdf';
    end
    % mean across subjects
    AccEns = cell(nSbj,1);
    AccROIs = cell(nSbj,1);
    for sbjitr = 1:nSbj
        AccEns{sbjitr} = res{sbjitr}.ensDec.profile_acc;
        AccROIs{sbjitr} = res{sbjitr}.mRoiDec.profile_acc_all(:,roiIdx);
    end
    ensacc = cellmean(AccEns);
    mroiacc = cellmean(AccROIs);
    
    % sort emotions based on ensemble
    [sorted,ord] = sort(ensacc,'descend');
    
    close all
    h = ffigure;
    for ix = 1:size(mroiacc,1)
        subplottight(r,c,o(ix),0.12);
        I = drawSurfaceMap(mroiacc(ord(ix),:),hcp360colMap,col);
        imagesc(I)
        set(gca,'FontSize',fsize)
        axis image off
        title(featname{ord(ix)})
    end
    
    suptitle(sprintf('Emotion score decoding accuracy: %s',scoreType));
    fname = [p.figdir,figname];
    pathtitle(fname)
    fprintf('%s\n',savprint(h,fname));
end
close all

%%
