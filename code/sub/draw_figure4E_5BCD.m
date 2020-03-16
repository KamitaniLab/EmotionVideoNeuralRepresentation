function draw_figure4E_5BCD(p,encRes)
%
% This code is for drawing figure 4E, 5BCD
%
%
%% settings
% use feature sets
compFeats = {'dimension','category','semantic','vision'};

load(p.roiInfPath,'roiInf');

nSbj = length(p.sbjID);

roiSets = {'individual','network17','subcortex'};

%% draw figures
for roisitr = 1:length(roiSets)
    switch roiSets{roisitr}
        case 'individual'
            roiDesc = roiInf.individual.roiDescrip;
            roiSet = roiInf.individual.roiName;
            yrang = [0,0.25];
            displayFeatType = {'VisSemCat'};
            fnames = {[p.figdir,'figure5B.pdf']};
        case 'network17'
            roiDesc = roiInf.network17.roiDescrip;
            roiSet = roiInf.network17.roiName;
            yrang = [0,0.2];
            displayFeatType = {'VisSemCat'};
            fnames = {[p.figdir,'figure5C.pdf']};
        case 'subcortex'
            roiDesc = roiInf.subcortex.roiDescrip;
            roiSet = roiInf.subcortex.roiName;
            yrang = [-0.02,0.06];
            displayFeatType = {'DimCat','VisSemCat'};
            fnames = {[p.figdir,'figure4E.pdf'],[p.figdir,'figure5D.pdf']};
    end
    nRois = length(roiSet);
    
    acc_all = cell(nSbj,1);
    for sbjitr = 1:nSbj
        %load metainf
        dpath = sprintf('%s%s/rois/%s_WholeBrain.mat',p.fmridir,p.sbjID{sbjitr},p.sbjID{sbjitr});
        load(dpath,'metainf');
        roiname = metainf.roiname;
        roiind_value = metainf.roiind_value;
        
        for scoritr = 1:length(compFeats)
            fset1 = compFeats{scoritr};
            acc1 = encRes.(fset1).pred_acc{sbjitr};
            
            % get values
            for ix = 1:nRois
                rois = sort(roiSet{ix});
                voxidx_cell = cell(length(rois),1);
                for ixx = 1:length(rois)
                    roiIdx = ismember(roiname,rois{ixx});
                    voxidx_cell{ixx} = any(roiind_value(roiIdx,:),1);
                end
                voxidx = any(merge(voxidx_cell),1);
                
                ac1 = acc1(voxidx);
                nanidx = isnan(ac1);
                nac1 = ac1(~nanidx);
                acc_all{sbjitr}{scoritr,ix} = nac1;
            end
        end
    end
    
    %% draw bars
    % display settings
    r = 6;
    c = 2;
    [r,c,o] = setrc2(r*c,'ltr',[r,c],[1,0]);
    fsize = 8;
    msize = 2;
    
    
    for dititr = 1:length(displayFeatType)
        
        close all
        h = ffigure;
        cnt = 0;
        
        switch displayFeatType{dititr}
            case 'DimCat'
                displayFeat = {'category','dimension'};
                cmap3('pr');
            case 'VisSemCat'
                displayFeat = {'vision','semantic','category'};
                cmap3('srb');
        end
        
        % assuming that compFeats = c:d, c:s,c:v
        acc_each_mu = cell(nSbj,1);
        acc_each_ci = cell(nSbj,1);
        scoreName = cell(length(compFeats),1);
        for scoritr = 1:length(compFeats)
            scoreName{scoritr} = compFeats{scoritr};
        end
        legName = displayFeat;
        for scoritr = 1:length(displayFeat)
            scIdx = ismember(scoreName,displayFeat{scoritr});
            for roitr = 1:length(roiSet)
                for sbjitr = 1:nSbj
                    [acc_each_ci{sbjitr}(scoritr,roitr),acc_each_mu{sbjitr}(scoritr,roitr)] = ciestim3(acc_all{sbjitr}{scIdx,roitr},2);
                end
            end
        end
        [acc_av_ci,acc_av_mu] = cellci(acc_each_mu);
        
        % accuracy
        for sbjitr = 1:nSbj
            cnt = cnt+1;
            subplottight(r,c,o(cnt),0.15);
            bar(acc_each_mu{sbjitr}','edgecolor','none')
            hold on
            shiftebar_h(acc_each_mu{sbjitr}',acc_each_ci{sbjitr}',length(displayFeat),'.k');
            set(gca,'FontSize',fsize)
            ylabel(sprintf('Mean acc (r)'))
            axname(roiDesc,1)
            xticklabel_rotate([],45)
            
            ylim(yrang)
            legend(legName,'Location','NorthEastOutside')
            hline(0,'-k')
            ffine(h)
            title(sprintf('Acc:%s',p.sbjID{sbjitr}),'FontSize',fsize)
            
        end
        % Average
        cnt = cnt+1;
        subplottight(r,c,o(cnt),0.15);
        hh = bar(acc_av_mu','edgecolor','none');
        hold on
        set(gca,'FontSize',fsize)
        hhh = get(hh,'Children');
        basePosi = zeros(length(displayFeat),length(roiSet));
        for scoritr = 1:length(displayFeat)
            basePosi(scoritr,:) = mean(get(hhh{scoritr},'Xdata'));
        end
        
        dif = 0;
        for sbjitr = 1:nSbj
            for roitr = 1:length(roiSet)
                xloc = repmat(basePosi(:,roitr)+dif,size(acc_each_mu{sbjitr}(:,roitr)',1),1);
                rand('seed',345)
                xloc = xloc+randn(size(xloc,1),1)/50;
                hhh = plot(xloc',acc_each_mu{sbjitr}(:,roitr),'-o','Color',[1,1,1]*0.4,'MarkerEdgeColor',[1,1,1]*0.2,'MarkerFaceColor',[1,1,1]*0.4,'MarkerSize',msize);
                set(gca,'FontSize',fsize)
            end
        end
        
        ylabel(sprintf('Mean acc (r)'))
        axname(roiDesc,1)
        xticklabel_rotate([],45)
        legend(legName,'Location','NorthEastOutside')
        hline(0,'-k')
        ffine(h)
        
        ylim(yrang)
        title(sprintf('%s','Average'),'FontSize',fsize)
        ffine(h)
        drawnow
        
        suptitle(sprintf('Mean encoding accuracy: %s, %s',roiSets{roisitr},displayFeatType{dititr}));
        pathtitle(fnames{dititr})
        fprintf('%s\n',savprint(h,fnames{dititr}));
    end
end


%%
close all