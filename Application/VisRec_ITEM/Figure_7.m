% ITEM Paper, Figure 7 (revised Figure 6)


clear
close all

%%% Step 1: load results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set subjects
load ../project_directories.mat
load ../subjects_all.mat

% set parameters
num_subj = numel(subj_ids);     % number of subjects
num_sess = 8;                   % number of sessions
num_sect = 48;                  % number of sectors

% preallocate correlations
RA = zeros(num_sess,num_sect,num_subj);
RS = zeros(num_sess,num_sect,num_subj);
RT = zeros(num_sess,num_sect,num_subj);

% for all subjects
for i = 1:num_subj
    
    % set directories
    subj_dir = strcat(bids_dir,'/','sub-',subj_ids{i});
    GLM_dir  = strcat(subj_dir,'/','anal','/','glms-item','/','glm-full');
    
    % load LS-A
    load(strcat(GLM_dir,'/','ITEM_dec_recon_LS_A','/','ITEM_sects-all_V1-all.mat'));
    for j = 1:num_sess, RA(j,:,i) = ITEM.Sess(j).CC; end;
    
    % load LS-S
    load(strcat(GLM_dir,'/','ITEM_dec_recon_LS_S','/','ITEM_sects-all_V1-all.mat'));
    for j = 1:num_sess, RS(j,:,i) = ITEM.Sess(j).CC; end;
    
    % load LS-T
    load(strcat(GLM_dir,'/','ITEM_dec_recon','/','ITEM_sects-all_V1-all.mat'));
    for j = 1:num_sess, RT(j,:,i) = ITEM.Sess(j).CC; end;
    
end;


%%% Step 2: analyze results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% compute averages
R_avg = zeros(3,num_sect,num_subj);
R_avg(1,:,:) = mean(RA,1);
R_avg(2,:,:) = mean(RS,1);
R_avg(3,:,:) = mean(RT,1);

% compute SEs
R_SEs = std(R_avg,[],3)./sqrt(num_subj);

% compute meta-averages
R_avg_avg = mean(R_avg,3);
R_avg_avg_avg = mean(R_avg_avg,2);

% test correlation coefficients
R_ft     = atan(R_avg_avg);
[h1, p1] = ttest(R_ft(3,:), R_ft(1,:));
[h2, p2] = ttest(R_ft(3,:), R_ft(2,:));
[h3, p3] = ttest(R_ft(2,:), R_ft(1,:));
fprintf('\n-> ITEM vs. LS-A: p = %e\n', p1);
fprintf('\n-> ITEM vs. LS-S: p = %e\n', p2);
fprintf('\n-> LS-S vs. LS-A: p = %e\n', p3);


%%% Step 3: display results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% display correlations
figure('Name', 'decoding accuracies', 'Color', [1 1 1], 'Position', [50 50 1200 1200]);
colors = 'rbg';
% colormap jet;
% cmap = colormap;
% cnum = size(cmap,1);

for h = 1:3
    subplot(3,1,h);
    hold on;
    bar(-1, 1, 'w');
    plot([-1, -1], [-1, +1], '--k', 'LineWidth', 2);
    plot([-1, -1], [-1, +1],  '-k', 'LineWidth', 2);
    for j = 1:num_sect
      % bar(j, R_avg_avg(j), 'FaceColor', cmap(1+round(((j-1)/(num_sect-1))*(cnum-1)),:));
        bar(j, R_avg_avg(h,j), colors(h));
        plot([(1-1), (num_sect+1)], [R_avg_avg_avg(h), R_avg_avg_avg(h)], '--k', 'LineWidth', 2);
        plot([j, j], [R_avg_avg(h,j)-R_SEs(h,j), R_avg_avg(h,j)+R_SEs(h,j)], '-k', 'LineWidth', 2);
        plot([j-1/4, j+1/4], [R_avg_avg(h,j)-R_SEs(h,j), R_avg_avg(h,j)-R_SEs(h,j)], '-k', 'LineWidth', 2);
        plot([j-1/4, j+1/4], [R_avg_avg(h,j)+R_SEs(h,j), R_avg_avg(h,j)+R_SEs(h,j)], '-k', 'LineWidth', 2);
    end;
    axis([(1-1), (num_sect+1), 0, 2/3]);
    set(gca, 'Box', 'On');
    set(gca, 'XTick', [1:num_sect]);
    set(gca, 'YTick', [0:0.1:0.6]);
    if h == 1, legend('mean over subjects', 'mean over sectors', 'SEM over subjects', 'Location', 'NorthEast'); end;
    if h == 3, xlabel('sector index', 'FontSize', 16); end;
    ylabel('correlation coefficient', 'FontSize', 16);
    if h == 1, title('LS-A', 'FontSize', 20); end;
    if h == 2, title('LS-S', 'FontSize', 20); end;
    if h == 3, title('ITEM', 'FontSize', 20); end;
end;