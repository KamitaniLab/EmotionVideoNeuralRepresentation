function draw_figure4D_S4D(p,encRes)
%
% This code is for drawing figure 4D and S4D
% Before drawing, slope estimation is conducted.
%
%% settings
compPairs = {
    {'category','dimension'}
    {'category','semantic'}
    {'category','vision'}
    };

load(p.roiInfPath,'roiInf');
p.roiDescrip = {roiInf.hcp360.roiDescrip{:},roiInf.individual.roiDescrip{:},roiInf.subcortex.roiDescrip{:}};
p.roiEnsemble = {roiInf.hcp360.roiDescrip{:},roiInf.subcortex.roiDescrip{:}}; % roi sets for ensemble decoders

nSbj = length(p.sbjID);

roiSets = {'WholeBrain','subcortex'};

hcpIdx = ~cellfun(@isempty,strfind(p.roiDescrip,'hcp180'));
subIdx = ismember(p.roiDescrip,{'thalamus','hyppocampus','hypothalamus','pallidum','brainstem','caudate','putamen','nuc_accum','amygdala','cerebellum'});
roiIdx = hcpIdx|subIdx;
nhcp = sum(hcpIdx);
nsub = sum(subIdx);

%% estimate slope values
clear Rslope
for roisitr = 1:length(roiSets)
    switch roiSets{roisitr}
        case 'WholeBrain'
            roiSet = {roiInf.hcp360.roiName{:},roiInf.subcortex.roiName{:}};
        case 'subcortex'
            roiSet = roiInf.subcortex.roiName;
    end
    nRois = length(roiSet);
    
    saveFname       = sprintf('%s/slope_estimate/%s/res.mat',p.savdir,roiSets{roisitr}); % res files
    setdir(fileparts(saveFname));
    
    if exist(saveFname,'file')
        %fprintf('Load:%s\n',saveFname)
        tmp =load(saveFname,'ac1_all','ac2_all','slopeDR','slopeDRp');
        Rslope.(roiSets{roisitr}) = tmp;
    else
        ac1_all = cell(nSbj,nRois);
        ac2_all = cell(nSbj,nRois);
        slopeDR = cell(length(compPairs),nSbj+1);
        slopeDRp = cell(length(compPairs),nSbj+1);
        for sbjitr = 1:nSbj
            %load metainf
            dpath = sprintf('%s%s/rois/%s_WholeBrain.mat',p.fmridir,p.sbjID{sbjitr},p.sbjID{sbjitr});
            load(dpath,'metainf');
            roiname = metainf.roiname;
            roiind_value = metainf.roiind_value;
            
            for ix = 1:nRois
                roiIdx = ismember(roiname,roiSet{ix});
                voxIdx = any(roiind_value(roiIdx,:),1);
                for compitr = 1:length(compPairs)
                    fset1 = compPairs{compitr}{1};
                    fset2 = compPairs{compitr}{2};
                    ac1 = encRes.(fset1).pred_acc{sbjitr}(:,voxIdx);
                    ac2 = encRes.(fset2).pred_acc{sbjitr}(:,voxIdx);
                    
                    nanidx = isnan(ac1+ac2);
                    nac1 = ac1(~nanidx);
                    nac2 = ac2(~nanidx);
                    ac1_all{sbjitr,ix}(:,compitr) = nac1;
                    ac2_all{sbjitr,ix}(:,compitr) = nac2;
                    [slopedr,int,st] = demingRegression(nac1',nac2',1,1,1,'both');
                    slopeDR{compitr,sbjitr}(ix) = slopedr;
                    slopeDRp{compitr,sbjitr}(ix) = st.pval;
                end
            end
        end
        
        % Pooled
        for compitr = 1:length(compPairs)
            for ix = 1:nRois
                ac1_pool = merge(ac1_all(:,ix),1);
                ac2_pool = merge(ac2_all(:,ix),1);
                [slopedr,int,st] = demingRegression(ac1_pool(:,compitr),ac2_pool(:,compitr),1,1,1,'both',100);
                slopeDR{compitr,nSbj+1}(ix) = slopedr;
                slopeDRp{compitr,nSbj+1}(ix) = st.pval;
            end
        end
        
        fprintf('Save:%s\n',saveFname)
        save(saveFname,'ac1_all','ac2_all','slopeDR','slopeDRp','-v7.3')
        
        Rslope.(roiSets{roisitr}).ac1_all = ac1_all;
        Rslope.(roiSets{roisitr}).ac2_all = ac2_all;
        Rslope.(roiSets{roisitr}).slopeDR = slopeDR;
        Rslope.(roiSets{roisitr}).slopeDRp = slopeDRp;
    end
end
clear tmp

%% draw figures

% display settings
r = 5;
c = 5;
[r,c,o] = setrc2(r*c,'ltr',[r,c],[1,1]);
fsize = 8;
msize = 4;%1.5;
msize_sub = 1.2;
msize_nsig = 1.2;
msize_cb = 1.6; % cerebellum
thres = 0; % deviation

axnam = {p.sbjID{:},'Pooled'};
nMultiComp = nhcp+nsub;

% get color index for hcp360 rois
tmp = load([p.rootPath,'code/util/hcp360colors.mat']);
hcpcols = tmp.vrank;

% get plot values
xd = cell(length(roiSets),length(compPairs));
yd = cell(length(roiSets),length(compPairs));
mp = cell(length(roiSets),length(compPairs));
for roisitr = 1:length(roiSets)
    slopeDR = Rslope.(roiSets{roisitr}).slopeDR;
    slopeDRp = Rslope.(roiSets{roisitr}).slopeDRp;
    
    for compitr = 1:length(compPairs)
        % angle
        mac = merge(slopeDR(compitr,:));
        mac = rad2deg(atan(mac))-45;
        % slope
        mp{roisitr,compitr} = merge(slopeDRp(compitr,:));
        
        hh = sinaplot(mac','MarkerType','.','edgecolor',[1,1,1]*0.5,'MarkerSize',msize,'mc',[],'ci',0);
        xd{roisitr,compitr} = get(hh,'XData');
        yd{roisitr,compitr} = get(hh,'YData');
        close all
    end
end

% draw figures
close all
h = ffigure;
cnt = 0;

for roisitr = 1:length(roiSets)
    switch roiSets{roisitr}
        case 'WholeBrain'
            roiSet = {roiInf.hcp360.roiName{:},roiInf.subcortex.roiName{:}};
        case 'subcortex'
            roiSet = roiInf.subcortex.roiName;
    end
    nRois = length(roiSet);
    
    for compitr = 1:length(compPairs)
        cnt = cnt + 1;
        subplottight(r,c,o(cnt),0.1);
        
        for i = 1:length(xd{roisitr,compitr})
            hold on
            % corrected by numbers of ROIs
            nonsig = (mp{roisitr,compitr}(i,:) >= (0.01/nMultiComp));
            sig = (mp{roisitr,compitr}(i,:) < (0.01/nMultiComp));
            catsig = (mp{roisitr,compitr}(i,:) < (0.01/nMultiComp)) & (yd{roisitr,compitr}{i} < thres);
            othersig = (mp{roisitr,compitr}(i,:) < (0.01/nMultiComp)) & (yd{roisitr,compitr}{i} > thres);
            
            sigIdx = find(sig);
            nonsigIdx = find(~sig);
            [s,ord] = sort(abs(yd{roisitr,compitr}{i}(sigIdx)),'ascend');
            sigIdx_sorted = sigIdx(ord);
            switch roiSets{roisitr}
                case {'WholeBrain'}
                    hcpcolssig = hcpcols(sig(1:nhcp));
                case {'subcortex'}
            end
            
            % draw non significant ROIs
            switch roiSets{roisitr}
                case {'WholeBrain'}
                    % cortex
                    hold on
                    plot(xd{roisitr,compitr}{i}(nonsigIdx(nonsigIdx<=nhcp)),yd{roisitr,compitr}{i}(nonsigIdx(nonsigIdx<=nhcp)),'o','MarkerEdgeColor',[1,1,1]*0.5,'MarkerSize',msize_nsig,'LineWidth',0.1);
                    % subcortex
                    if any(nonsigIdx(nonsigIdx>nhcp))
                        subnonsigIdx = nonsigIdx(nonsigIdx>nhcp);
                        for ixxx = 1:length(subnonsigIdx)
                            plot(xd{roisitr,compitr}{i}(subnonsigIdx(ixxx)),yd{roisitr,compitr}{i}(subnonsigIdx(ixxx)),plotMarker2(subnonsigIdx(ixxx)-nhcp),'MarkerEdgeColor',[1,1,1]*0.5,'MarkerSize',msize_nsig,'LineWidth',0.1);
                        end
                    end
                case {'subcortex'}
                    % cortex
                    hold on
                    % subcortex
                    subnonsigIdx = nonsigIdx;
                    for ixxx = 1:length(subnonsigIdx)
                        plot(xd{roisitr,compitr}{i}(subnonsigIdx(ixxx)),yd{roisitr,compitr}{i}(subnonsigIdx(ixxx)),plotMarker2(subnonsigIdx(ixxx)),'MarkerEdgeColor',[1,1,1]*0.5,'MarkerSize',msize_nsig*2,'LineWidth',0.1);
                    end
            end
            % draw significant ROIs
            for ix = 1:length(sigIdx_sorted)
                % cortex
                switch roiSets{roisitr}
                    case {'WholeBrain'}
                        if sigIdx_sorted(ix) <= nhcp
                            plot(xd{roisitr,compitr}{i}(sigIdx_sorted(ix)),yd{roisitr,compitr}{i}(sigIdx_sorted(ix)),'.','MarkerEdgeColor',plotColor(hcpcolssig(ord(ix)),nRois/2,'Spectral_r'),'MarkerSize',msize,'LineWidth',0.1);
                        elseif (sigIdx_sorted(ix) == nhcp+8) || (sigIdx_sorted(ix) == nhcp+9) % +x
                            % subcortex
                            plot(xd{roisitr,compitr}{i}(sigIdx_sorted(ix)),yd{roisitr,compitr}{i}(sigIdx_sorted(ix)),plotMarker2(sigIdx_sorted(ix)-nhcp),...
                                'MarkerEdgeColor',[1,1,1]*0.2,...
                                'MarkerFaceColor',[1,1,1]*0.2,...
                                'MarkerSize',msize_cb,'LineWidth',0.2);
                        else
                            % subcortex
                            plot(xd{roisitr,compitr}{i}(sigIdx_sorted(ix)),yd{roisitr,compitr}{i}(sigIdx_sorted(ix)),plotMarker2(sigIdx_sorted(ix)-nhcp),...
                                'MarkerEdgeColor',[1,1,1]*0.2,...
                                'MarkerFaceColor',[1,1,1]*0.2,...
                                'MarkerSize',msize_sub,'LineWidth',0.1);
                        end
                        
                    case {'subcortex'}
                        if (sigIdx_sorted(ix) == 8)||(sigIdx_sorted(ix) == 9) % cerebellum
                            % subcortex
                            plot(xd{roisitr,compitr}{i}(sigIdx_sorted(ix)),yd{roisitr,compitr}{i}(sigIdx_sorted(ix)),plotMarker2(sigIdx_sorted(ix)),...
                                'MarkerEdgeColor',[1,1,1]*0.2,...
                                'MarkerFaceColor',[1,1,1]*0.2,...
                                'MarkerSize',msize_cb*2,'LineWidth',1);
                        else
                            % subcortex
                            plot(xd{roisitr,compitr}{i}(sigIdx_sorted(ix)),yd{roisitr,compitr}{i}(sigIdx_sorted(ix)),plotMarker2(sigIdx_sorted(ix)),...
                                'MarkerEdgeColor',[1,1,1]*0.2,...
                                'MarkerFaceColor',[1,1,1]*0.2,...
                                'MarkerSize',msize_sub*2,'LineWidth',0.1);
                        end
                        
                end
            end
            text(i,20,sprintf('%d',sum(othersig)),'FontSize',fsize)
            text(i,-20,sprintf('%d',sum(catsig)),'FontSize',fsize)
        end
        set(gca,'FontSize',fsize)
        hline(thres,'-k')
        ylabel(sprintf('Deviation from the parity'))
        axname(axnam,1)
        axname([-30:15:30],2,[-30:15:30])
        xticklabel_rotate([],45)
        for ix = (length(p.sbjID)+1)
            text(ix,25,compPairs{compitr}{2}(1:3),'FontSize',fsize)
            text(ix,-25,compPairs{compitr}{1}(1:3),'FontSize',fsize)
        end
        ylim([-30,30])
        title(sprintf('%s: %s vs %s:',roiSets{roisitr},compPairs{compitr}{1},compPairs{compitr}{2}))
        ffine(h)
    end
    
end
suptitle(sprintf('Distribution of deviations from the parity'));
fname = [p.figdir,'figure4D_S4D.pdf'];
pathtitle(fname)
fprintf('%s\n',savprint(h,fname));


%%
close all
