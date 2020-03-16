function encRes = summary_encoding_results(p)
%
% - This code is written for summarizing whole brain encoding analysis.
%
% called from emotion2020_analysis_BATCH
%

%% load results
for sbjitr = 1:length(p.sbjID)
    sbjID = p.sbjID{sbjitr};
    
    for scoritr = 1:length(p.scoreTypes)
        scoreType = p.scoreTypes{scoritr};
        % set suffix
        suffix_full = [sbjID,'/',p.suffix_enc];
        dataFname            = sprintf('%s/%s/%s/res.mat',p.savdir,scoreType,suffix_full); % res files
        fprintf('Load:%s\n',dataFname)
        tmp = load(dataFname,'r_best');
        
        encRes.(scoreType).pred_acc{sbjitr} = tmp.r_best.pred_acc;
        encRes.(scoreType).pred_acc_cv{sbjitr} = tmp.r_best.pred_acc_cv;
        encRes.(scoreType).iden_acc{sbjitr} = tmp.r_best.iden_acc;
        clear tmp
    end
end

%%
