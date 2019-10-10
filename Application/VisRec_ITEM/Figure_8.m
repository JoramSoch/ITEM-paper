% ITEM Paper, Figure 8


clear
% close all

%%% Step 1: load results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set subjects
load ../project_directories.mat
load ../subjects_all.mat

% set parameters
num_subj = numel(subj_ids);     % number of subjects
num_sess = 8;                   % number of sessions
num_sect = 48;                  % number of sectors
num_trls = 100;                 % number of trials

% preallocate correlations
R = zeros(num_sess,num_sect,num_subj);

% for all subjects
for i = 1:num_subj
    
    % set directories
    subj_dir = strcat(bids_dir,'/','sub-',subj_ids{i});
    GLM_dir  = strcat(subj_dir,'/','anal','/','glms-item','/','glm-full');
    
    % load correlations
    load(strcat(GLM_dir,'/','ITEM_dec_recon','/','ITEM_sects-all_V1-all.mat'));
    for j = 1:num_sess, R(j,:,i) = ITEM.Sess(j).CC; end;
    
end;


%%% Step 2: analyze results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute max, min, med
r_max = max(max(max(R)));
r_min = min(min(min(abs(R))));
r_med = median([R(:); 1]);

% locate max, min, med
[j_max, k_max, i_max] = ind2sub([num_sess, num_sect, num_subj], find(R==r_max));
[j_min, k_min, i_min] = ind2sub([num_sess, num_sect, num_subj], find(abs(R)==r_min));
[j_med, k_med, i_med] = ind2sub([num_sess, num_sect, num_subj], find(R==r_med));

% load actual and reconstructed values (max)
subj_dir = strcat(bids_dir,'/','sub-',subj_ids{i_max});
GLM_dir  = strcat(subj_dir,'/','anal','/','glms-item','/','glm-full');
load(strcat(GLM_dir,'/','ITEM_dec_recon','/','ITEM_sects-all_V1-all.mat'));
yt_max = ITEM.Sess(j_max).Yt(:,k_max);
yr_max = ITEM.Sess(j_max).Yr(:,k_max);

% load actual and reconstructed values (min)
subj_dir = strcat(bids_dir,'/','sub-',subj_ids{i_min});
GLM_dir  = strcat(subj_dir,'/','anal','/','glms-item','/','glm-full');
load(strcat(GLM_dir,'/','ITEM_dec_recon','/','ITEM_sects-all_V1-all.mat'));
yt_min = ITEM.Sess(j_min).Yt(:,k_min);
yr_min = ITEM.Sess(j_min).Yr(:,k_min);

% load actual and reconstructed values (med)
subj_dir = strcat(bids_dir,'/','sub-',subj_ids{i_med});
GLM_dir  = strcat(subj_dir,'/','anal','/','glms-item','/','glm-full');
load(strcat(GLM_dir,'/','ITEM_dec_recon','/','ITEM_sects-all_V1-all.mat'));
yt_med = ITEM.Sess(j_med).Yt(:,k_med);
yr_med = ITEM.Sess(j_med).Yr(:,k_med);


%%% Step 3: display results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% display correlations
figure('Name', 'exemplary reconstructions', 'Color', [1 1 1], 'Position', [50 50 1200 1200]);

for h = 1:3
    subplot(3,1,h);
    hold on;
    if h == 1                   % min
        plot([1:num_trls], yt_min, '-or', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
        plot([1:num_trls], yr_min, '-ob', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
    end;
    if h == 2                   % med
        plot([1:num_trls], yt_med, '-or', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
        plot([1:num_trls], yr_med, '-ob', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
    end;
    if h == 3                   % max
        plot([1:num_trls], yt_max, '-or', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
        plot([1:num_trls], yr_max, '-ob', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'w');
    end;    
    xlim([(1-1), (num_trls+1)]);
    ylim([-3.25, +3.25]);
    set(gca, 'Box', 'On');
    if h == 1, legend('presented contrast', 'reconstructed contrast', 'Location', 'NorthEast'); end;
    if h == 3, xlabel('trial index', 'FontSize', 16); end;
    ylabel('contrast value', 'FontSize', 16);
    if h == 1
        title('worst reconstruction', 'FontSize', 20);
        text(num_trls/2, -3, sprintf('r = %0.2f (Subject %d, Session %d, Sector %d)', r_min, i_min, j_min, k_min), 'FontSize', 12, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom');
    end;
    if h == 2
        title('medium reconstruction', 'FontSize', 20);
        text(num_trls/2, -3, sprintf('r = %0.2f (Subject %d, Session %d, Sector %d)', r_med, i_med, j_med, k_med), 'FontSize', 12, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom');
    end;
    if h == 3,
        title('best reconstruction', 'FontSize', 20);
        text(num_trls/2, -3, sprintf('r = %0.2f (Subject %d, Session %d, Sector %d)', r_max, i_max, j_max, k_max), 'FontSize', 12, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom');
    end;
end;