function draw_figureS1AB(p,decRes)
%
% This code is for drawing figure S1AB
%
%
%% settings
% score index
dimIdx = 1:14;
catIdx = 15:48;
emoIdx = 1:48;

% get roi Idx
indivIdx = ismember(p.roiDescrip,{'VC','TPJ','IPL','PC','STS','TE','MTC','Insula','DLPFC','DMPFC','VMPFC','ACC','OFC'});
subIdx = ismember(p.roiDescrip,{'thalamus','hyppocampus','hypothalamus','pallidum','brainstem','caudate','putamen','nuc_accum','amygdala','cerebellum'});

nIndivRois = sum(indivIdx);
nSubRois = sum(subIdx);

nSbj = length(p.sbjID);

% get feature names
load([p.featdir,'category.mat']);
catname = L.featname';
load([p.featdir,'dimension.mat']);
dimname = L.featname';


roiSetTypes = {'individual','sub'};

%% draw figure

% display settings
c = 2;
r = 6;
[r,c,o] = setrc2(r*c,'ltr',[r,c],[0]);
fsize = 5;
msize = 2;
col = 'Spectral';

AccIndiv = cell(nSbj,1);
AccSub = cell(nSbj,1);
EnsCatAcc = cell(nSbj,1);
EnsDimAcc = cell(nSbj,1);
for sbjitr = 1:nSbj
    catRes = decRes.category{sbjitr};
    dimRes = decRes.dimension{sbjitr};
    
    AccIndiv{sbjitr} = [catRes.mRoiDec.profile_acc_all(:,indivIdx);dimRes.mRoiDec.profile_acc_all(:,indivIdx)];
    AccSub{sbjitr} = [catRes.mRoiDec.profile_acc_all(:,subIdx);dimRes.mRoiDec.profile_acc_all(:,subIdx)];
    EnsCatAcc{sbjitr} = catRes.ensDec.profile_acc;
    EnsDimAcc{sbjitr} = dimRes.ensDec.profile_acc;
end

% draw figures
close all
h = ffigure;
cnt = 0;

% draw individual ROI accuracy
for roisitr = 1:length(roiSetTypes)
    switch roiSetTypes{roisitr}
        case 'individual'
            acc = AccIndiv;
            roinames = {'VC','TPJ','IPL','PC','STS','TE','MTC','Insula','DLPFC','DMPFC','VMPFC','ACC','OFC'};
            nrois = nIndivRois;
            pos = ([0.6308,0.6923,0.7538,0.8154,0.8769,0.9385,1.0000,1.0615,1.1231,1.1846,1.2462,1.3077,1.3692]-1)/3;
        case 'sub'
            acc = AccSub;
            roinames = {'thalamus','hyppocampus','hypothalamus','pallidum','brainstem','caudate','putamen','nuc_accum','amygdala','cerebellum'};
            nrois = nSubRois;
            pos = ([0.6400,0.7200,0.8000,0.8800,0.9600,1.0400,1.1200,1.2000,1.2800,1.3600]-1)/3;
    end
    
    % individual subject accuracy
    for sbjitr = 1:nSbj
        [s,ord_dim] = sort(EnsDimAcc{sbjitr},'descend');
        [s,ord_cat] = sort(EnsCatAcc{sbjitr},'descend');
        cnt = cnt+1;
        subplottight(r,c,o(cnt),0.1);
        for ix = 1:nrois
            plot((1:length(emoIdx))+pos(ix),acc{sbjitr}([ord_cat;ord_dim+34],ix),plotMarker2(ix),...
                'MarkerFaceColor','none',...
                'MarkerEdgeColor',plotColor(ix,nrois,col),...
                'MarkerSize',msize)
            hold on
            switch roiSetTypes{roisitr}
                case 'individual'
                    if ix == nrois
                        plot((1:length(emoIdx))+pos(ix),acc{sbjitr}([ord_cat;ord_dim+34],ix),plotMarker2(ix+1),...
                            'MarkerFaceColor','none',...
                            'MarkerEdgeColor',plotColor(ix,nrois,col),...
                            'MarkerSize',msize)
                    end
            end
        end
        ylim([-0.1,0.6])
        set(gca,'FontSize',fsize)
        hh = hline(0.1:0.1:0.7,'-k');
        set(hh,'Color',[1,1,1]*0.2)
        hline(0,'-k');
        vline(length(catIdx)+0.5,'-k');
        hline(length(catIdx)+0.5,'-k');
        ylabel('Corr.Coeff.')
        axname(strrep([catname(ord_cat),dimname(ord_dim)],'_',' '),1)
        xticklabel_rotate([],45)
        legend(roinames,'Location','NorthEastOutSide')
        title(sprintf('%s:%s',roiSetTypes{roisitr},p.sbjID{sbjitr}))
        ffine(h)
    end
    
    % average accuracy
    [acc_ci,acc_mu] = cellci(acc);
    cnt = cnt+1;
    subplottight(r,c,o(cnt),0.1);
    [s,ord_dim] = sort(cellmean(EnsDimAcc),'descend');
    [s,ord_cat] = sort(cellmean(EnsCatAcc),'descend');
    for ix = 1:nrois
        plot((1:length(emoIdx))+pos(ix),acc_mu([ord_cat;ord_dim+34],ix),plotMarker2(ix),...
            'MarkerFaceColor','none',...
            'MarkerEdgeColor',plotColor(ix,nrois,col),...
            'MarkerSize',msize)
        switch roiSetTypes{roisitr}
            case 'individual'
                if ix == nrois
                    plot((1:length(emoIdx))+pos(ix),acc_mu([ord_cat;ord_dim+34],ix),plotMarker2(ix+1),...
                        'MarkerFaceColor','none',...
                        'MarkerEdgeColor',plotColor(ix,nrois,col),...
                        'MarkerSize',msize)
                end
        end
        hold on
    end
    ylim([-0.1,0.6])
    set(gca,'FontSize',fsize)
    hh = hline(0.1:0.1:0.7,'-k');
    set(hh,'Color',[1,1,1]*0.2)
    hline(0,'-k');
    vline(length(catIdx)+0.5,'-k');
    hline(length(catIdx)+0.5,'-k');
    ylabel('Corr.Coeff.')
    axname(strrep([catname(ord_cat),dimname(ord_dim)],'_',' '),1)
    xticklabel_rotate([],45)
    legend(roinames,'Location','NorthEastOutSide')
    title(sprintf('%s:%s',roiSetTypes{roisitr},'Average'))
    ffine(h)
end

suptitle(sprintf('Decoding accuracy for emotion scores from individual ROIs'));
fname = [p.figdir,'figureS1AB.pdf'];
pathtitle(fname)
fprintf('%s\n',savprint(h,fname));

%%
close all