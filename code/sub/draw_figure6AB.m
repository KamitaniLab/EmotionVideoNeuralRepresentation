function draw_figure6AB(p)
%
% This code is for drawing figure 6AB
% Umap analysis (umap_on_brain_activity_distance.ipynb) must be finished before executing this script.
%
%
%% settings
mapPath = [p.rootPath,'res/umap_on_brain_activity_distance/'];

% color preparation
load([p.featdir,'category.mat']);
ck27from34 = [1:11,13,15:17,19:20,22:25,27:32]; % select 27 dimensions as Cowen & Keltner (2017)
cat27 = L.feat(:,ck27from34);
cat27(p.dupidx,:) = [];
ckcol = merge(cmap4('ck27',0),1);

c = cat27;
c(c<0) = 0;
cols = (c*ckcol);
cols = cols/max(cols(:));
cols = cols*2;
cols(cols(:)>1) = 1;
cols = cols/max(cols(:));
cols(cols(:)>=0.9999) = 0.9999;

% sort sample presentation order
nscore = normal(cat27,2);
nscore = nscore./repmat(sum(nscore,2),1,27);
H = zeros(2181,1);
for i = 1:2181
    H(i) = entropy0(nscore(i,:));
end
[sorted,ord] = sort(H,'descend');


%% Umap parameters
n_neighbour = 20;
min_dist = 0.5;
seed = 1;
spread = 10;

%% draw figure
% figure settings
fsize = 5;
r = 1;
c = 8;
[r,c,o] = setrc2(r*c,'ltd',[r,c],[0,1]);
msize = 1;


close all
h = ffigure;
cnt = 0;

sbjTypes = {p.sbjID{:},'Average'};

for sbjitr = 1:length(sbjTypes)
    sbjType = sbjTypes{sbjitr};
    
    paramName    = sprintf('nn%d_md%.1f_sp%d_seed%d_%s',n_neighbour,min_dist,spread,seed,sbjType);
    dataFname    = sprintf('%s%s.mat',mapPath,paramName);
    
    tmp = load(dataFname);
    
    cnt = cnt + 1;
    subplottight(r,c,o(cnt),0.1);
    for i = 1:size(tmp.dist_umap,1)
        plot(tmp.dist_umap(ord(i),1),tmp.dist_umap(ord(i),2),'o','MarkerFaceColor',cols(ord(i),:),'MarkerEdgeColor',cols(ord(i),:),'MarkerSize',msize)
        hold on
    end
    axis image off
    title(sprintf('%s:%s',strrep(paramName,'_','-'),sbjType),'FontSize',fsize)
end
suptitle(sprintf('Umap on brain activity pattern distance with 27 category color'));
fname = [p.figdir,'figure6A.pdf'];
pathtitle(fname)
fprintf('%s\n',savprint(h,fname));

close all




%% draw figure 6B

% color settings
dim = load([p.featdir,'dimension.mat']);
dim14 = dim.L.feat;
dim14(p.dupidx,:) = [];

dim28c = load([p.featdir,'dim28continuous.mat']);
dim28c.L.feat(p.dupidx,:) = [];
pdim28 = dim28c.L.feat(:,1:14);
ndim28 = dim28c.L.feat(:,15:28);
pvad = dim28c.L.feat(:,[14,2,7]);
nvad = dim28c.L.feat(:,[14,2,7]+14);
[sorted,ord] = sort(pdim28(:,end)+ndim28(:,end),'ascend'); % valence sort

sbjTypes = {p.sbjID{:},'Average'};

% figure settings
close all
h = ffigure;
cnt = 0;
fsize = 5;
r = 4;
c = 8;
[r,c,o] = setrc2(r*c,'ltd',[r,c],[1]);
msize = 1;

for sbjitr = 1:length(sbjTypes)
    sbjType = sbjTypes{sbjitr};
    scores = dim28c.L.feat;
    scname = dim28c.L.featname;
    
    
    paramName    = sprintf('nn%d_md%.1f_sp%d_seed%d_%s',n_neighbour,min_dist,spread,seed,sbjType);
    dataFname    = sprintf('%s%s.mat',mapPath,paramName);
    tmp = load(dataFname);
    
    % rgb valence arousal dominance
    rgbcol = eye(3);
    pvadcols = (pvad*rgbcol);
    pvadcols = pvadcols/max(pvadcols(:));
    scalefactor = 1.3;
    pvadcols = pvadcols*scalefactor;
    pvadcols(pvadcols(:)>1) = 1;
    pvadcols = pvadcols/max(pvadcols(:));
    pvadcols(pvadcols(:)>=0.9999) = 0.9999;
    
    nvadcols = (nvad*rgbcol);
    nvadcols = nvadcols/max(nvadcols(:));
    scalefactor = 1.3;
    nvadcols = nvadcols*scalefactor;
    nvadcols(nvadcols(:)>1) = 1;
    nvadcols = nvadcols/max(nvadcols(:));
    nvadcols(nvadcols(:)>=0.9999) = 0.9999;
    
    cnt = cnt + 1;
    subplottight(r,c,o(cnt),0.1);
    for i = 1:size(tmp.dist_umap,1)
        plot(tmp.dist_umap(ord(i),1),tmp.dist_umap(ord(i),2),'o','MarkerFaceColor',pvadcols(ord(i),:),...
            'MarkerEdgeColor',pvadcols(ord(i),:),'MarkerSize',msize)
        hold on
    end
    axis image off
    title(sprintf('%s','vad_high'),'FontSize',fsize)
    
    cnt = cnt + 1;
    subplottight(r,c,o(cnt),0.1);
    for i = 1:size(tmp.dist_umap,1)
        plot(tmp.dist_umap(ord(i),1),tmp.dist_umap(ord(i),2),'o','MarkerFaceColor',nvadcols(ord(i),:),...
            'MarkerEdgeColor',nvadcols(ord(i),:),'MarkerSize',msize)
        hold on
    end
    axis image off
    title(sprintf('%s','vad_low'),'FontSize',fsize)
    
end

suptitle(sprintf('Umap on brain activity pattern distance with 3 dimension color'));
fname = [p.figdir,'figure6B.pdf'];
pathtitle(fname)
fprintf('%s\n',savprint(h,fname));


%% draw distributions of individual scores

% color settings
cat = load([p.featdir,'category.mat']);
cat27 = cat.L.feat(:,[1:11,13,15:17,19:20,22:25,27:32]);
cat27(p.dupidx,:) = [];

dim28c = load([p.featdir,'dim28continuous.mat']);
dim28 = dim28c.L.feat;
dim28(p.dupidx,:) = [];

sbjTypes = {'Average'};
colTypes = {'cat','dim'};

% sort scores
for sbjitr = 1:length(sbjTypes)
    sbjType = sbjTypes{sbjitr};
    for colitr = 1:length(colTypes)
        colType = colTypes{colitr};
        
        switch colType
            case 'cat'
                scores = cat27;
                scname = cat.L.featname([1:11,13,15:17,19:20,22:25,27:32]);
            case 'dim'
                scores = dim28;
                scname = dim28c.L.featname;
        end
        [sorted,ord] = sort(scores,'ascend');
        nScores = size(scores,2);
        seq = 1:nScores;
        
        % figure settings
        close all
        h = ffigure;
        cnt = 0;
        fsize = 5;
        r = 6;
        c = 9;
        [r,c,o] = setrc2(r*c,'ltr',[r,c],[1]);
        msize = 1;
        
        paramName    = sprintf('nn%d_md%.1f_sp%d_seed%d_%s',n_neighbour,min_dist,spread,seed,sbjType);
        dataFname    = sprintf('%s%s.mat',mapPath,paramName);
        tmp = load(dataFname);
        
        for ix = 1:nScores
            
            % empty irrelevant colors
            switch colType
                case 'cat'
                    ckcol = merge(cmap4('ck27',0),1);
                    cs = cat27;
                    scalefactor = 2;
                case 'dim'
                    ckcol = merge(cmap4('ck28c',0),1);
                    cs = dim28;
                    scalefactor = 1.3;
            end
            cs(:,seq ~= ix) = 0;
            cs(cs<0) = 0;
            cols = (cs*ckcol);
            cols = cols/max(cols(:));
            cols = cols*scalefactor;
            cols(cols(:)>1) = 1;
            cols = cols/max(cols(:));
            cols(cols(:)>=0.9999) = 0.9999;
            
            cnt = cnt + 1;
            subplottight(r,c,o(cnt),0.1);
            for i = 1:size(tmp.dist_umap,1)
                plot(tmp.dist_umap(ord(i,ix),1),tmp.dist_umap(ord(i,ix),2),'o','MarkerFaceColor',cols(ord(i,ix),:),'MarkerEdgeColor',cols(ord(i,ix),:),'MarkerSize',msize)
                hold on
            end
            axis image off
            title(sprintf('%s',scname{ix}),'FontSize',fsize)
        end
        
        suptitle(sprintf('Umap on brain activity pattern distance with individual score color:%s:%s',sbjType,colType));
        fname = [p.figdir,'figure6AB_',sbjType,colType,'.pdf'];
        pathtitle(fname)
        fprintf('%s\n',savprint(h,fname));
        
    end
    
end


%%
close all
