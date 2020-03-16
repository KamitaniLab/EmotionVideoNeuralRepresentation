function rsaRes = summary_rsa_results(p)
%
% - This code is written for summarizing representational similarity analysis.
%
% called from emotion2020_analysis_BATCH
%

%% load results
for sbjitr = 1:length(p.sbjID)
    sbjID = p.sbjID{sbjitr};
    
    for scoritr = 1:length(p.scoreTypes)
        scoreType = p.scoreTypes{scoritr};
        switch scoreType
            case {'dimension','category','semantic','vision'};
            otherwise
                continue
        end
        fprintf('Load RS results:%s:%s\n',sbjID,scoreType)
        
        for roitr = 1:length(p.roiDescrip)
            roi = p.roiDescrip{roitr};
            suffix_full = [sbjID,'/',p.suffix_rsa,'/',roi];
            dataFname            = sprintf('%s/%s/%s/res.mat',p.savdir,scoreType,suffix_full); % res files
            tmp = load(dataFname,'res');
            
            rsaRes.(scoreType).rs{sbjitr}(roitr) = tmp.res.corr;
        end
        clear tmp
    end
end

%%
