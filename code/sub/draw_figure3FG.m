function draw_figure3FG(p,decRes)
%
% - This code is written for performing emotion identification analysis
%   based on region-wise decoding accuracy
%
%
% called from emotion2020_analysis_BATCH
%

%% settings
% score index
dimIdx = 1:14;
catIdx = 15:48;
emoIdx = 1:48;

colcat = [216,31,73]/255;
coldim = [250,120,210]/255;

nSbj = length(p.sbjID);

% identification parameters
nrepetitions = 100;

%% draw figures

% display settings
c = 5;
r = 4;
[r,c,o] = setrc2(r*c,'ltr',[r,c],[1,1]);
msize = 2;
fsize = 8;

dataFname = sprintf('%s/emotion_identification/result_%d.mat',p.savdir,nrepetitions);
if ~exist(dataFname,'file')
    fprintf('The results were not yet prepared. Please run perform_emotion_identification_analysis.m\n')
else
    load(dataFname);

    % prepare within- across- accuracy map similarity
    npairs = nSbj*(nSbj-1);
    within_cat = cell(npairs,1);
    within_dim = cell(npairs,1);
    across_cat = cell(npairs,1);
    across_dim = cell(npairs,1);
    for npitr = 1:npairs
        within_cat{npitr} = diag(CatSimMat_all{npitr});
        within_dim{npitr} = diag(DimSimMat_all{npitr});
        across_cat{npitr} = mean(rmvDiag(CatSimMat_all{npitr}),2);
        across_dim{npitr} = mean(rmvDiag(DimSimMat_all{npitr}),2);
    end
    
    figure;
    mtc = {cellmean(within_cat),cellmean(across_cat),cellmean(within_dim),cellmean(across_dim)};
    hh = sinaplot(mtc);
    xd = get(hh,'XData');
    yd = get(hh,'YData');
    close
    
    % draw figure
    close all
    h = ffigure;
    cnt = 0;
    
    % draw within- across- accuracy map similarity
    cnt = cnt +1;
    subplottight(r,c,o(cnt),0.1);
    mt = [[mean(cellmean(within_cat),1),mean(cellmean(across_cat),1)];[mean(cellmean(within_dim),1),mean(cellmean(across_dim),1)]];
    bar(mt);
    hold on
    
    cols = cell(4,1);
    cols{1} = colcat;
    cols{2} = colcat;
    cols{3} = coldim;
    cols{4} = coldim;
    diffs = [-1,-2,-3,-4];
    adds = [0.8572,1.1429,1.8572,2.1429];
    
    plot([[(xd{1}+diffs(1))/4+adds(1)];[(xd{2}+diffs(2))/4+adds(2)]],[yd{1};yd{2}],'-','Color',[1,1,1]/2,'LineWidth',0.5)
    hold on
    plot([[(xd{3}+diffs(3))/4+adds(3)];[(xd{4}+diffs(4))/4+adds(4)]],[yd{3};yd{4}],'-','Color',[1,1,1]/2,'LineWidth',0.5)
    hold on
    
    for jj = 1:length(xd)
        plot((xd{jj}+diffs(jj))/4+adds(jj),yd{jj},'o','MarkerFaceColor',cols{jj},'MarkerEdgeColor',[1,1,1]*0.5,'MarkerSize',msize);
        hold on
    end
    set(gca,'FontSize',fsize)
    axname({'Category','Dimension'})
    ffine(h)
    axis square
    ylabel('Correlation coefficient')
    title('Within/Across Corr.')
    
    % identification accuracy
    cnt = cnt+1;
    subplottight(r,c,o(cnt),0.1);
    dif = -0.1;
    hh = errorbar_h((1:(length(catIdx)-1))'+dif,CatIdenAcc_mu_all,CatIdenAcc_ci_all);
    set(gca,'FontSize',fsize)
    set(hh,'Color',colcat);
    AnnotationOff(hh,1:2:length(hh));
    hold on
    dif = 0.1;
    hh = errorbar_h((1:(length(dimIdx)-1))+dif,DimIdenAcc_mu_all,DimIdenAcc_ci_all);
    set(hh,'Color',coldim);
    AnnotationOff(hh,1:2:length(hh));
    plot(100./(2:34),'--k')
    xlim([0,35])
    ylim([0,100])
    ylabel('Identification accuracy (%)')
    xlabel('# of candidates')
    legend({'Category','Dimension.','Chance'},'Location','NorthEast')
    axname(2:4:34,1,1:4:33)
    ffine(h)
    axis square
    title('Emotion identificaion accuracy')
    
    suptitle(sprintf('Region-wise accuracy map similarity, and emotion identification performance'));
    fname = [p.figdir,'figure3FG.pdf'];
    pathtitle(fname)
    fprintf('%s\n',savprint(h,fname));
    
end

%%
