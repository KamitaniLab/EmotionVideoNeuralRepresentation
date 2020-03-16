function draw_figureS5(p,decRes,rsaRes)
%
% This code is for drawing figure S5
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

analysisTypes = {'decSingle','decDouble','rsaSingle','rsaDouble'};

%% draw figure S5
% figure settings
fsize = 8;
r = 6;
c = 6;
[r,c,o] = setrc2(r*c,'ltr',[r,c],[1]);

% initialize
close all
h = ffigure;
cnt = 0;

for analitr = 1:length(analysisTypes)
    
    switch analysisTypes{analitr}
        case 'decSingle'
            col = 'hot';
            range = [50,90];
            scoreTypes = {'category','dimension','vision','semantic'};
            AccROIs = cell(nSbj,length(scoreTypes));
            for scoritr = 1:length(scoreTypes)
                for sbjitr = 1:nSbj
                    AccROIs{sbjitr,scoritr} = mean(decRes.(scoreTypes{scoritr}){sbjitr}.mRoiDec.iden_acc_all(:,roiIdx));
                end
            end
            
        case 'rsaSingle'
            col = 'hot';
            range = [0,0.1];
            scoreTypes = {'category','dimension','vision','semantic'};
            AccROIs = cell(nSbj,length(scoreTypes));
            for scoritr = 1:length(scoreTypes)
                for sbjitr = 1:nSbj
                    AccROIs{sbjitr,scoritr} = rsaRes.(scoreTypes{scoritr}).rs{sbjitr}(roiIdx);
                end
            end
            
        case 'decDouble'
            col = 'RdBu_r';
            range = [-5,5];
            scoreTypes = {
                {'category','vision'}
                {'category','semantic'}
                {'dimension','vision'}
                {'dimension','semantic'}
                };
            AccROIs = cell(nSbj,length(scoreTypes));
            for scoritr = 1:length(scoreTypes)
                for sbjitr = 1:nSbj
                    AccROIs{sbjitr,scoritr} = mean(decRes.(scoreTypes{scoritr}{1}){sbjitr}.mRoiDec.iden_acc_all(:,roiIdx))-...
                        mean(decRes.(scoreTypes{scoritr}{2}){sbjitr}.mRoiDec.iden_acc_all(:,roiIdx));
                end
            end
        case 'rsaDouble'
            col = 'RdBu_r';
            range = [-0.03,0.03];
            scoreTypes = {
                {'category','vision'}
                {'category','semantic'}
                {'dimension','vision'}
                {'dimension','semantic'}
                };
            AccROIs = cell(nSbj,length(scoreTypes));
            for scoritr = 1:length(scoreTypes)
                for sbjitr = 1:nSbj
                    AccROIs{sbjitr,scoritr} = rsaRes.(scoreTypes{scoritr}{1}).rs{sbjitr}(roiIdx)-...
                        rsaRes.(scoreTypes{scoritr}{2}).rs{sbjitr}(roiIdx);
                end
            end
            
    end
    
    for scoritr = 1:length(scoreTypes)
        cnt = cnt+1;
        mroiacc = cellmean(AccROIs(:,scoritr));
        
        subplottight(r,c,o(cnt),0.12);
        I = drawSurfaceMap(mroiacc,hcp360colMap,col,range);
        imagesc(I)
        set(gca,'FontSize',fsize)
        axis image off
        switch analysisTypes{analitr}
            case {'decSingle','rsaSingle'}
                title(sprintf('%s:%s',analysisTypes{analitr},scoreTypes{scoritr}))
            case {'decDouble','rsaDouble'}
                title(sprintf('%s:%s-%s',analysisTypes{analitr},scoreTypes{scoritr}{1},scoreTypes{scoritr}{2}))
        end
    end
    
end

suptitle(sprintf('Decoding identification and RSA performances'));
fname = [p.figdir,'figureS5.pdf'];
pathtitle(fname)
fprintf('%s\n',savprint(h,fname));

close all

%%
