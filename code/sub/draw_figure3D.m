function draw_figure3D(p,decRes)
%
% This code is for drawing figure 3ABC and S1C
% Umap analysis (umap_on_emotion_scores.ipynb) must be finished before executing this script.
%
%
%% settings
mapPath = [p.rootPath,'res/umap_on_emotion_scores/'];

% prepare label inf. for score based coloring
dpath = sprintf('%s%s/rois/%s_WholeBrain.mat',p.fmridir,p.sbjID{1},p.sbjID{1});
load(dpath,'metainf');
label_index = metainf.Label;
inds_all = 1:length(label_index);
inds_all = inds_all(~ismember(label_index,p.dupidx));
load([p.featdir,'category.mat']);

% color preparation
ck27from34 = [1:11,13,15:17,19:20,22:25,27:32]; % select 27 dimensions as Cowen & Keltner (2017)
cat27 = L.feat(label_index(inds_all),ck27from34);
ckcol = merge(cmap4('ck27',0),1);
c = cat27-repmat(mean(cat27,1),size(cat27,1),1);
c(c<0) = 0;
cols = (c*ckcol);
cols = cols/max(cols(:));
cols = cols*2.2;
cols(cols(:)>1) = 1;
cols = cols/max(cols(:));
cols(cols(:)>=0.9999) = 0.9999;


nSbj = length(p.sbjID);

%% Umap parameters
n_neighbour = 10;
min_dist = 1.0;
seed = 1;
spread = 3;

%% figure 3D

% figure settings
fsize = 5;
r = 1;
c = nSbj+3;
[r,c,o] = setrc2(r*c,'ltd',[r,c],[0,1]);
msize = 1;

close all
h = ffigure;
cnt = 0;

% decoded
x = zeros(nSbj,1);
y = zeros(nSbj,1);
mu = zeros(nSbj,1);
for sbjitr = 1:nSbj
    sbj = p.sbjID{sbjitr};
    paramName    = sprintf('nn%d_md%.1f_sp%d_seed%d_%s',n_neighbour,min_dist,spread,seed,sbj);
    dataFname    = sprintf('%s%s.mat',mapPath,paramName);
    
    tmp = load(dataFname);
    ed = diag(sqdist(tmp.dec_umap',tmp.cat_umap'));
    [sorted,ord] = sort(ed,'descend');
    
    cnt = cnt + 1;
    subplottight(r,c,o(cnt),0.1);
    for i = 1:size(tmp.dec_umap,1)
        plot(tmp.dec_umap(ord(i),1),tmp.dec_umap(ord(i),2),'o','MarkerFaceColor',cols(ord(i),:),'MarkerEdgeColor',cols(ord(i),:),'MarkerSize',msize)
        hold on
    end
    axis image off
    x(sbjitr) = fcorr(tmp.dec_umap(:,1),tmp.cat_umap(:,1));
    y(sbjitr) = fcorr(tmp.dec_umap(:,2),tmp.cat_umap(:,2));
    mu(sbjitr) = (x(sbjitr)+y(sbjitr))/2;
    title(sprintf('%s\n[x,y,mu]=[%.3f,%.3f,%.3f]',strrep(paramName,'_','-'),x(sbjitr),y(sbjitr),mu(sbjitr)),'FontSize',fsize)
end

% true
% every subject data contain true score data
paramName    = sprintf('nn%d_md%.1f_sp%d_seed%d',n_neighbour,min_dist,spread,seed);
dataFname    = sprintf('%s%s_%s.mat',mapPath,paramName,p.sbjID{1});
tmp = load(dataFname);

cnt = cnt + 1;
subplottight(r,c,o(cnt),0.1);
for i = 1:size(tmp.cat_umap,1)
    plot(tmp.cat_umap(i,1),tmp.cat_umap(i,2),'o','MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:),'MarkerSize',msize)
    hold on
end
axis image off

% show averaged accuracy
xa = mean(x,1);
ya = mean(y,1);
mua = mean(mu,1);
title(sprintf('True:%sAverage\n:[x,y,mu]=[%.3f,%.3f,%.3f]',strrep(paramName(1:end-3),'_','-'),xa,ya,mua),'FontSize',fsize)

suptitle(sprintf('Umap on true and decoded feature'));
fname = [p.figdir,'figure3D.pdf'];
pathtitle(fname)
fprintf('%s\n',savprint(h,fname));


close all

%% plot with individual scores for checking emotion positioning

% figure settings
fsize = 5;
r = 6;
c = 9;
[r,c,o] = setrc2(r*c,'ltr',[r,c],[1]);
msize = 1;

close all
h = ffigure;
cnt = 0;
paramName    = sprintf('nn%d_md%.1f_sp%d_seed%d_%s',n_neighbour,min_dist,spread,seed,p.sbjID{1});
dataFname    = sprintf('%s%s.mat',mapPath,paramName);
tmp = load(dataFname);

[sorted,ord] = sort(cat27,'ascend');
nScores = size(cat27,2);
scname = L.featname([1:11,13,15:17,19:20,22:25,27:32]);
seq = 1:nScores;

for ix = 1:nScores
    
    % ignore irrelevant colors
    cs = cat27-repmat(mean(cat27,1),size(cat27,1),1);
    cs(:,seq ~= ix) = 0;
    cs(cs<0) = 0;
    cols = (cs*ckcol);
    cols = cols/max(cols(:));
    cols = cols*2.2;
    cols(cols(:)>1) = 1;
    cols = cols/max(cols(:));
    cols(cols(:)>=0.9999) = 0.9999;
    
    cnt = cnt + 1;
    subplottight(r,c,o(cnt),0.1);
    for i = 1:size(tmp.cat_umap,1)
        plot(tmp.cat_umap(ord(i,ix),1),tmp.cat_umap(ord(i,ix),2),'o','MarkerFaceColor',cols(ord(i,ix),:),'MarkerEdgeColor',cols(ord(i,ix),:),'MarkerSize',msize)
        hold on
    end
    axis image off
    
    title(sprintf('True:%s\n',scname{ix}),'FontSize',fsize)
end

suptitle(sprintf('Umap :%s',strrep(paramName(1:end-3),'_','-')));
fname = [p.figdir,'figure3D_supp.pdf'];
pathtitle(fname)
fprintf('%s\n',savprint(h,fname));


%%
close all
