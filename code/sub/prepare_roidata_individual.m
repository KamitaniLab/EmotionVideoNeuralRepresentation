function p = prepare_roidata_individual(param)
%
% This script was written for creating preprocessed MRI data for cowen &
% keltner movie clip presentation data.
%
% written by Tomoyasu Horikawa horikawa-t@atr.jp 2020/03/11
%
%
% Called from prepare_roidata.m
% 
% Required functions:
%  setdir.m, merge.m
% 
% Required libraries:
% BrainDecoderToolbox2
% 

%% get parameters
p = param;
p.mfilename = mfilename;

%% get path to individual data
% set path to data
fmri_mat_file = cell(p.nSession,1);
for i = 1:p.nSession
    fmri_mat_file{i} = sprintf('%sfmri_%s_Session%d.h5',p.fmridir_org,p.sbj,i);
end

% load individual dataset
dataSet_all  = cell(p.nSession,1);
metaData_all = cell(p.nSession,1);
for i = 1:p.nSession
    fprintf('Load:Session%d\n',i)
    [dataSet_all{i},metaData_all{i}] = load_data(fmri_mat_file{i});
end

%% integrate multiple data and save individual ROI data separately
% set meta info.
metaInfList = {'Session','Run','Block','Label'};

for ritr = 1:length(p.roiDescrip)
    
    % set path & log information
    p.rois{ritr,1} = p.roiDescrip{ritr};
    saveFname = [p.fmridir_rois,p.sbj '_'  p.rois{ritr,1},'.mat'];
    saveFnameChk = [p.logdir_sbj,p.sbj '_'  p.rois{ritr,1},'_log.txt'];
    p.savdir_sbj_rois{ritr,1} = saveFname;
    p.savdir_sbj_roisChk{ritr,1} = saveFnameChk;
    
    if ~exist(saveFnameChk,'file')
        fprintf('Start:%s\n',saveFnameChk)
        
        % save log file
        tmp = [];
        save(saveFnameChk,'tmp', '-ascii')
        %delete(saveFnameChk)
        
        % extract data
        braindat_all = cell(p.nSession,1);
        clear metainf_all metainf
        for ditr = 1:p.nSession
            
            % get voxel data
            bdataKey = 'VoxelData';
            vdind = ismember(metaData_all{ditr}.key,bdataKey);
            
            % get xyz coordinate
            if ditr == 1
                xind = ismember(metaData_all{ditr}.key,'voxel_x');
                yind = ismember(metaData_all{ditr}.key,'voxel_y');
                zind = ismember(metaData_all{ditr}.key,'voxel_z');
                x = metaData_all{ditr}.value(xind,metaData_all{ditr}.value(vdind,:) == 1);
                y = metaData_all{ditr}.value(yind,metaData_all{ditr}.value(vdind,:) == 1);
                z = metaData_all{ditr}.value(zind,metaData_all{ditr}.value(vdind,:) == 1);
                metainf_all.xyz = [x;y;z];
            end
            
            % get meta info
            for metitr = 1:length(metaInfList)
                mind = ismember(metaData_all{ditr}.key,metaInfList{metitr});
                metainf_all.(metaInfList{metitr}){ditr} = dataSet_all{ditr}(:,metaData_all{ditr}.value(mind,:) == 1);
                
                % increment for mutiple dataset
                switch metaInfList{metitr}
                    case {'Run','Block','Session'}
                        if ditr == 1
                            cnt.(metaInfList{metitr}) = 0;
                        else
                            cnt.(metaInfList{metitr}) = metainf_all.(metaInfList{metitr}){ditr-1}(end);
                        end
                        metainf_all.(metaInfList{metitr}){ditr} = metainf_all.(metaInfList{metitr}){ditr} + cnt.(metaInfList{metitr});
                    case 'Label'
                        labelIdx = 3; % stimulus label index
                        metainf_all.(metaInfList{metitr}){ditr} = metainf_all.(metaInfList{metitr}){ditr}(:,labelIdx);
                end
            end
            
            % get ROI (only conduct for 1st dataset) assuming the same data size across sessions.
            % specify ROI descriptions
            if ditr == 1
                roiind_desc = ~cellfun(@isempty,strfind(metaData_all{ditr}.description,'ROI'));
                roiind_value = metaData_all{ditr}.value(roiind_desc,metaData_all{ditr}.value(vdind,:) == 1) == 1;
                roiind_key = metaData_all{ditr}.key(roiind_desc);
                roiskey = strrep(strrep(roiind_key,'_r_rh','_rh'),'_r_lh','_lh');

                % use rois
                useroiind_cell = cell(length(p.roiName{ritr}),1);
                rois = strrep(strrep(p.roiName{ritr},'_r_rh','_rh'),'_r_lh','_lh');
                for roinamitr = 1:length(p.roiName{ritr})
                    useroiind_cell{roinamitr} = strcmp(roiskey,rois{roinamitr});
                end
                useroiind_vec = any(merge(useroiind_cell,2),2);
                useroiind_num = find(useroiind_vec);
                for i = useroiind_num
                    fprintf('use: %s\n',roiind_key{i})
                end
                fprintf('\n')
                
                % exclude rois
                exroiind_cell = cell(length(p.exRoiName{ritr}),1);
                rois = strrep(strrep(p.exRoiName{ritr},'_r_rh','_rh'),'_r_lh','_lh');
                for roinamitr = 1:length(p.exRoiName{ritr})
                    exroiind_cell{roinamitr} = strcmp(roiskey,rois{roinamitr});
                end
                exroiind_vec = any(merge(exroiind_cell,2),2);
                exroiind_num = find(exroiind_vec);
                for i = exroiind_num
                    fprintf('exclude: %s\n',roiind_key{i})
                end
                fprintf('\n')
                
                % preserve voxel index for individual ROIs
                if ~isempty(exroiind_vec)
                    voxind_all = roiind_value(useroiind_vec & ~exroiind_vec,:);
                else
                    voxind_all = roiind_value(useroiind_vec,:);
                end
                voxind_use = any(voxind_all,1);
                
                % extract used voxel coordinate
                metainf_all.xyz = metainf_all.xyz(:,voxind_use);
                
                % keep ROI-voxel info
                metainf_all.roiname = roiind_key;
                metainf_all.roiind_value = roiind_value(:,voxind_use);
                metainf_all.roiname_used = p.roiName{ritr};
                metainf_all.roiname_excluded = p.exRoiName{ritr};
                %metainf_all.voxind_all = voxind_all(:,voxind_use);
            end
            
            % extract ROIs
            braindat_tmp = dataSet_all{ditr}(:,metaData_all{ditr}.value(vdind,:) == 1);
            braindat_all{ditr} = braindat_tmp(:,voxind_use);
            clear braindat_tmp
        end
        
        % merge all
        braindat = single(merge(braindat_all,1));
        metainf.xyz = metainf_all.xyz;
        metainf.roiname = metainf_all.roiname;
        metainf.roiind_value = metainf_all.roiind_value == 1;
        metainf.roiname_used = metainf_all.roiname_used';
        metainf.roiname_excluded = metainf_all.roiname_excluded;
        %metainf.voxind_all = metainf_all.voxind_all == 1;
        for metitr = 1:length(metaInfList)
            metainf.(metaInfList{metitr}) = merge(metainf_all.(metaInfList{metitr}),1);
        end
        metainf.label_type = {'stimID'};
        metainf
        
        % save preprocessed data
        fprintf('Save:%s\n',saveFname)
        save(saveFname,'braindat','metainf','-v7.3')
        clear braindat_all
    else
        fprintf('Skip:%s\n',saveFnameChk)
    end
end

end % end function

%%


