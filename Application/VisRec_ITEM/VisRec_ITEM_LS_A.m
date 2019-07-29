% Visual Reconstruction: complete data analysis (***least squares, all***)
% _
% This script performs data analysis for visual reconstruction paradigm.
% 
% Author: Joram Soch
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 30/11/2018, 15:55 / 03/12/2018, 17:40 /
%         05/12/2018, 17:20 / 19/12/2018, 23:30 /
%         19/02/2019, 11:35 / 01/04/2019, 11:30


%%% Step 0: Set up analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

% steps to perform
thx2do = [1,5];
sub2do = [1:4];

% load directories
cd('..');
load project_directories.mat

% set subjects
if numel(sub2do) == 1, subj_file = 'subjects_1st.mat'; end;
if numel(sub2do) == 3, subj_file = 'subjects_2-4.mat'; end;
if numel(sub2do) == 4, subj_file = 'subjects_all.mat'; end;

% set models
MS_file = 'model_spaces/glms-item.mat';

% set ROIs
ROI_imgs = {'V1-left',  'V1-right', ...
            'V2a-left', 'V2b-left', 'V2a-right', 'V2b-right', ...
            'V3a-left', 'V3b-left', 'V3a-right', 'V3b-right', ...
            'V4-left',  'V4-right'}';
ROIs_VR  = [1 2];

% set estimation mode, regressors of no interest, modulator transformation
est_mod = 'DCT';                % 'KXY'  = filter design and data
                                % 'DCT'  = use discrete cosine set
reg_noi = 'incl';               % 'incl' = include into 2nd-lvl model
                                % 'excl' = exclude from 2nd-lvl model
PM_trans = 'lin';               % 'lin'  = linear     : x-0.5
                                % 'exp'  = exponential: 10^(x-1) & mc
                                % 'log'  = logarithmic: log10((10-10^0.1)*x+10^0.1)) & mc
                                %     (x = 0, 1/3, 2/3, 1)
                                %    (mc = mean-centering)

% load subjects and models
load(subj_file);
load(MS_file);
num_subj = numel(subj_ids);
num_ROIs = numel(ROI_imgs);

% set other parameters
num_sess = 8;
num_scan = 220;
num_trls = 100;
num_sect = 48;

% set model space
ind_sect = [1:num_sect];
num_mods = numel(ind_sect);


%%% Step 1: Estimate first-level model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(1,thx2do)

cd(strcat(tool_dir,'VisRec_ITEM/LS-A_LS-S/'));
    
% display message
fprintf('\n\n-> Estimate first-level model:\n');

% for all subjects
for i = 1:num_subj
    
    % load SPM.mat
    fprintf('   - Subject %s (%d out of %d) ... ', subj_ids{i}, i, num_subj);
    SPM_mat = strcat(bids_dir,'/','sub-',subj_ids{i},'/anal/','glms-',MS_name,'/','glm-',GLM_names{1},'/','SPM.mat');
    load(SPM_mat);
    
    % estimate model
    ITEM_est_1st_lvl_LS_A(SPM, est_mod, [2 3]);
    fprintf('successful!\n');

end;

cd(tool_dir);

end;


%%% Step 5: Estimate second-level model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(5,thx2do)

cd(strcat(tool_dir,'VisRec_ITEM/LS-A_LS-S/'));
    
% display message
fprintf('\n\n-> Estimate second-level model:\n');
    
% for all subjects
for i = 1:num_subj
    
    % load SPM
    fprintf('   - Subject %s (%d out of %d): ', subj_ids{i}, i, num_subj);
    SPM_mat = strcat(bids_dir,'/','sub-',subj_ids{i},'/anal/','glms-',MS_name,'/','glm-',GLM_names{1},'/','SPM.mat');
    load(SPM_mat);
    
    % configure decoding
    c   = [0, ones(1,num_sect)];
    con = 'sects-all';
    
    % configure region
    ROI_nii = strcat(SPM.swd,'/','ITEM_est_1st_lvl','/','ROI_','V1-all','_LBF_all-m0.nii');
    reg     = 'V1-all';
    
    % perform decoding
    ITEM_dec_recon_LS_A(SPM, ROI_nii, c, con, reg);
    fprintf('done.\n');
    
end;

cd(tool_dir);

end;


%%% Step 6: Display visual reconstruction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(6,thx2do)

% display message
fprintf('\n\n-> Display visual reconstruction:\n');

% preallocate correlations
r = zeros(num_sess,num_sect,num_subj);

% load correlations
for i = 1:num_subj
    fprintf('   - Subject %s (%d out of %d): ', subj_ids{i}, i, num_subj);
    con = 'sects-all';
    reg = 'V1-all';
    ITEM_mat = strcat(bids_dir,'/','sub-',subj_ids{i},'/anal/','glms-',MS_name,'/','glm-',GLM_names{1},'/','ITEM_dec_recon_LS_A','/','ITEM_',con,'_',reg,'.mat');
    load(ITEM_mat);
    r(:,:,i) = vertcat(ITEM.Sess.CC);
    fprintf('done.\n');
end;

% open figure
figure('Name','Visual Reconstruction (LS-A)','Color',[1 1 1],'Position',[50 50 1600 900])

% display correlations
for i = 1:num_subj
    CC = r(:,:,i);
    subplot(1,num_subj,i);
    bar([1:num_sess], mean(CC,2), 'b');
    axis([(1-0.5), (num_sess+0.5), 0, 0.5]);
    set(gca,'Box','On');
    set(gca,'XTick',[1:num_sess]);
    xlabel('recording session', 'FontSize', 12);
    ylabel('average correlation', 'FontSize', 12);
    title(sprintf('Subject %d: %s', i, subj_ids{i}), 'FontSize', 12);
end;

end;

cd(strcat(tool_dir,'VisRec_ITEM/'));