% Visual Reconstruction: searchlight-based analysis
% _
% This script performs visual reconstruction using searchlight-based ITEM.
% 
% Author: Joram Soch
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 10/05/2019, 15:20


%%% Step 0: Set up analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

% steps to perform
thx2do = [5];
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

% load subjects and models
load(subj_file);
load(MS_file);
num_subj = numel(subj_ids);

% set other parameters
num_sess = 8;
num_scan = 220;
num_trls = 100;
num_sect = 48;

% set model space
ind_sect = [1:48];
num_mods = numel(ind_sect);


%%% Step 5: Estimate second-level model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismember(5,thx2do)

% display message
fprintf('\n\n-> Estimate second-level model:\n');
    
% for all subjects
for i = 1:num_subj
    
    % load SPM
    SPM_mat = strcat(bids_dir,'/','sub-',subj_ids{i},'/anal/','glms-',MS_name,'/','glm-',GLM_names{1},'/','SPM.mat');
    load(SPM_mat);
    
    % configure decoding
    rad = 9;
    c   = [0, ones(1,num_sect)];
    con = 'sects-all';
    
    % for each searchlight
    fprintf('   - Subject %s (%d out of %d):\n', subj_ids{i}, i, num_subj);
    fprintf('     - searchlight-based ITEM: ');
    ITEM_dec_recon_SL(SPM, rad, c, con)
    fprintf('done.\n');
    
end;

end;