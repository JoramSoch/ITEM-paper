% Visual Reconstruction: complete data analysis
% _
% This script performs data analysis for visual reconstruction paradigm.
% 
% Author: Joram Soch
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 30/11/2018, 15:55 / 03/12/2018, 17:40 /
%         05/12/2018, 17:20 / 19/12/2018, 23:30


%%% Step 0: Set up analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

% steps to perform
thx2do = [0:6];
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
ind_sect = [1:48];
num_mods = numel(ind_sect);


%%% Step 0: Estimate standard model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(0,thx2do)

% prepare model estimation
modelling_meta_batch(subj_file, MS_file, 'prepare')

% perform model estimation
modelling_meta_batch(subj_file, MS_file, 'perform')

end;


%%% Step 1: Estimate first-level model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(1,thx2do)

% display message
fprintf('\n\n-> Estimate first-level model:\n');

% for all subjects
for i = 1:num_subj
    
    % load SPM.mat
    fprintf('   - Subject %s (%d out of %d) ... ', subj_ids{i}, i, num_subj);
    SPM_mat = strcat(bids_dir,'/','sub-',subj_ids{i},'/anal/','glms-',MS_name,'/','glm-',GLM_names{1},'/','SPM.mat');
    load(SPM_mat);
    
    % estimate model
    ITEM_est_1st_lvl(SPM, est_mod, [2 3]);
    fprintf('successful!\n');

end;

end;


%%% Step 2: Extract regions of interest %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(2,thx2do)

% display message
fprintf('\n\n-> Load regions of interest:\n');
    
% for all subjects
for i = 1:num_subj

    % load SPM.mat
    SPM_mat = strcat(bids_dir,'/','sub-',subj_ids{i},'/anal/','glms-',MS_name,'/','glm-',GLM_names{1},'/','SPM.mat');
    load(SPM_mat);
    
    % load mask.nii
    M = MA_load_mask(SPM);
    R = struct([]);

    % preallocate ROIs
    m_ind   = [];
    num_vox = zeros(num_ROIs,1);
    
    % load all ROIs
    fprintf('   - Subject %s (%d out of %d): ROIs ', subj_ids{i}, i, num_subj);
    for j = 1:num_ROIs
        % load ROI image
        fprintf('%d, ', j);
        roi_str = strcat(bids_dir,'/','sub-',subj_ids{i},'/anat/','rsub-',subj_ids{i},'_roi-',ROI_imgs{j},'.nii');
        roi_hdr = spm_vol(roi_str);
        roi_img = spm_read_vols(roi_hdr);
        roi_img = reshape(roi_img,[1 prod(roi_hdr.dim)]);
        roi_ind = find(roi_img~=0 & M~=0); % all ROI voxels within mask
        % get ROI voxels
        num_vox(j) = numel(roi_ind);
        R(j).Y_ind = [(numel(m_ind)+1):(numel(m_ind)+num_vox(j))];
        R(j).m_ind = roi_ind;
        R(j).v     = num_vox(j);
        m_ind = [m_ind, roi_ind];
    end;
    clear roi_str roi_hdr roi_img roi_ind
    fprintf('end.\n');
    
    % save ROI data
    filename = strcat(SPM.swd,'/','ITEM_est_1st_lvl','/','ROIs.mat');
    save(filename, 'R');
    
end;
    
end;


%%% Step 3: Bayesian model assessment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(3,thx2do)

cd(strcat(tool_dir,'VisRec_ITEM/'));

% for all subjects
for i = 1:num_subj
    
    % load SPM, GLM, ROIs
    SPM_mat = strcat(bids_dir,'/','sub-',subj_ids{i},'/anal/','glms-',MS_name,'/','glm-',GLM_names{1},'/','SPM.mat');
    load(SPM_mat);
    GLM_mat = strcat(SPM.swd,'/','ITEM_est_1st_lvl','/','GLM1.mat');
    load(GLM_mat);
    ROI_mat = strcat(GLM1.swd,'/','ROIs.mat');
    load(ROI_mat);
    
    % get ROI voxels
    m_ind = horzcat(R(ROIs_VR).m_ind);
    GLM_Y = struct([]);
    GLM_X = struct([]);
    
    % load gamma estimates
    v = numel(m_ind);
    G = ITEM_load_gammas(SPM, m_ind);
    
    % perform ReML analysis
    % for j = 1:num_sess
    %     X0 = GLM1.Sess(j).T(:,[1, (1+num_sect+1):end]);
    %     [GLM_Y(j).V, GLM_Y(j).s2] = ITEM_GLM_ReML(G(GLM1.Sess(j).t,:), X0, eye(GLM1.tr(j)), GLM1.Sess(j).U);
    %     [GLM_Y(j).V, GLM_Y(j).s2] = ITEM_GLM_ReML(G(GLM1.Sess(j).t,:), GLM1.Sess(j).T, eye(GLM1.tr(j)), GLM1.Sess(j).U, ...
    %                                               sprintf('VisRec_ITEM, Step 3: ReML estimation for session %d', j));
    % end;
    
    % assemble data and models
    for j = 1:num_sess
        GLM_Y(j).Y = G(GLM1.Sess(j).t,:);
        GLM_Y(j).P = inv(GLM1.Sess(j).U);
        for k = 0:num_mods
            if k == 0
                GLM_X(num_sect+1,j).X0 = GLM1.Sess(j).T(:,[1, (1+num_sect+1):end]);
            else
                GLM_X(k,j).X = GLM1.Sess(j).T(:,[1, 1+ind_sect(k), (1+num_sect+1):end]);
            end;
        end;
    end;
    
    % preallocate results
    GLM    = GLM_Y;
    cvLME0 = zeros(1,numel(m_ind));
    cvLMEs = zeros(num_mods,numel(m_ind));
    fprintf('\n\n-> Calculate cvLMEs for Subject %s (%d out of %d):\n', subj_ids{i}, i, num_subj);
    
    % for all models
    for k = 0:num_mods
        % null model
        if k == 0
            fprintf('   - Estimating model %d: null model ... ', k);
            for j = 1:num_sess, GLM(j).X = GLM_X(num_sect+1,j).X0; end;
            [cvLME0, oosLME] = ME_GLM_cvLME(GLM);
        end;
        % other models
        if k > 0
            fprintf('   - Estimating model %d: sector %d ... ', k, ind_sect(k));
            for j = 1:num_sess, GLM(j).X = GLM_X(k,j).X; end;
            [cvLMEs(k,:), oosLME] = ME_GLM_cvLME(GLM);
        end;
        fprintf('successful!\n');
    end;
    
    % save results
    clear oosLME
    filename = strcat(GLM1.swd,'/','LMEs.mat');
    save(filename, 'GLM_Y', 'GLM_X', 'cvLMEs', 'cvLME0');
    
end;

end;


%%% Step 4: Bayesian model comparison %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(4,thx2do)

% display message
fprintf('\n\n-> Calculate cvLBFs for all vs. m0:\n');

% for all subjects
for i = 1:num_subj
    
    % load SPM & GLM
    SPM_mat = strcat(bids_dir,'/','sub-',subj_ids{i},'/anal/','glms-',MS_name,'/','glm-',GLM_names{1},'/','SPM.mat');
    GLM_mat = strcat(bids_dir,'/','sub-',subj_ids{i},'/anal/','glms-',MS_name,'/','glm-',GLM_names{1},'/','ITEM_est_1st_lvl','/','GLM1.mat');
    load(SPM_mat);
    load(GLM_mat);
    
    % load ROIs & LMEs
    ROI_mat = strcat(GLM1.swd,'/','ROIs.mat');
    LME_mat = strcat(GLM1.swd,'/','LMEs.mat');
    load(ROI_mat);
    load(LME_mat);
    
    % load mask
    M = MA_load_mask(SPM);
    H = MA_init_header(SPM, false);
    
    % compare all models against null model
    cvLFE_all  = ME_MF_LFE(cvLMEs);
    LBF_all_m0 = cvLFE_all - cvLME0;
    
    % for each ROI
    fprintf('   - Subject %s (%d out of %d): ROIs ', subj_ids{i}, i, num_subj);
    roi_all = [];
    for j = ROIs_VR
        % determine voxels with most evidence in ROI
        fprintf('%d, ', j);
        roi_mat = [R(j).m_ind', LBF_all_m0(R(j).Y_ind)'];
        roi_mat = flipud(sortrows(roi_mat, 2));
        roi_ind = roi_mat(1:num_sect, 1);
        roi_val = roi_mat(1:num_sect, 2);
        roi_all = [roi_all; roi_ind, roi_val];
        % write restricted ROI image to disk
        roi_img = NaN(size(M));
        roi_img(roi_ind) = roi_val;
        H.fname   = strcat(GLM1.swd,'/','ROI_',ROI_imgs{j},'_LBF_all-m0.nii');
        H.descrip = sprintf('VisRec_ITEM: voxels with most evidence for visual processing in %s', ROI_imgs{j});
        spm_write_vol(H, reshape(roi_img, H.dim));
    end;
    clear roi_mat roi_ind roi_val roi_img
    
    % for all ROIs
    fprintf('all, ');
    roi_img = NaN(size(M));
    roi_img(roi_all(:,1)) = roi_all(:,2);
    H.fname   = strcat(GLM1.swd,'/','ROI_','V1-all','_LBF_all-m0.nii');
    H.descrip = sprintf('VisRec_ITEM: voxels with most evidence for visual processing in %s', 'V1-all');
    spm_write_vol(H, reshape(roi_img, H.dim));
    fprintf('done.\n');
    
end;

end;


%%% Step 5: Estimate second-level model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(5,thx2do)

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
    ITEM_dec_recon(SPM, ROI_nii, c, con, reg);
    fprintf('done.\n');
    
end;

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
    ITEM_mat = strcat(bids_dir,'/','sub-',subj_ids{i},'/anal/','glms-',MS_name,'/','glm-',GLM_names{1},'/','ITEM_dec_recon','/','ITEM_',con,'_',reg,'.mat');
    load(ITEM_mat);
    r(:,:,i) = vertcat(ITEM.Sess.CC);
    fprintf('done.\n');
end;

% open figure
figure('Name','Visual Reconstruction','Color',[1 1 1],'Position',[50 50 1600 900])

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