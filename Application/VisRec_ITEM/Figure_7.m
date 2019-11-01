% ITEM Paper, Figure 7


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

% compute first averages
R_avg_avg = mean(R_avg,3);
R_avg_SEs = std(R_avg,[],3)./sqrt(num_subj);

% compute second averages
R_avg_avg_avg = mean(R_avg_avg,2);
R_avg_avg_SEs = std(R_avg_avg,[],2)./sqrt(num_sect);

% test correlation coefficients
R_ft     = atan(R_avg_avg);
[h1, p1] = ttest(R_ft(3,:), R_ft(1,:));
[h2, p2] = ttest(R_ft(3,:), R_ft(2,:));
[h3, p3] = ttest(R_ft(2,:), R_ft(1,:));
fprintf('\n-> ITEM vs. LS-A: p = %e', p1);
fprintf('\n-> ITEM vs. LS-S: p = %e', p2);
fprintf('\n-> LS-S vs. LS-A: p = %e\n\n', p3);


%%% Step 3: display results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd('StimDiscLib/');

% display correlations
figure('Name', 'decoding accuracies', 'Color', [1 1 1], 'Position', [50 50 1200 1200]);
colors = 'rbg';
method = {'LS-A', 'LS-S', 'ITEM'};
bx     = [1, 3; 2, 3; 1, 2];
by     = [0.45, 0.40, 0.35];
dx     =  0.1;
dy     =  0.025;

for h = 1:3
    % correlation means
    subplot(3,3,(h-1)*3+1);
    plot_dartboard('-k', 2);
    plot_dataset(R_avg_avg(h,:), 'jet', 256, [0, 2/3], 'None', '-k', 2);
    text((-1-0.25), 0, method{h}, 'FontSize', 20, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle', 'Rotation', 90);
    text((+1+0.45), 0, '±', 'FontSize', 28, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Baseline');
    if h == 1, title('mean', 'FontSize', 20); end;
    % correlation SEMs
    subplot(3,3,(h-1)*3+2);
    plot_dartboard('-k', 2);
    plot_dataset(R_avg_SEs(h,:), 'jet', 256, [0, 1/6], 'None', '-k', 2);
    if h == 1, title('SEM', 'FontSize', 20); end;
    % sector indices
    if h == 2
        subplot(3,3,(h-1)*3+3);
        plot_dartboard('-k', 2);
        plot_dataset(0.2*ones(1,48), 'gray', 256, [-1, +1], 'None', '-k', 2);
        plot_sectors('w', 8, 'bold');
        title('sectors', 'FontSize', 20); 
    end;
    % average averages
    if h == 3
        subplot(3,3,(h-1)*3+3);
        hold on;
        for i = 1:3
            bar(i, R_avg_avg_avg(i), colors(i));
            plot([i, i], [R_avg_avg_avg(i)-R_avg_avg_SEs(i), R_avg_avg_avg(i)+R_avg_avg_SEs(i)], '-k', 'LineWidth', 2);
            plot([i-dx, i+dx], [R_avg_avg_avg(i)-R_avg_avg_SEs(i), R_avg_avg_avg(i)-R_avg_avg_SEs(i)], '-k', 'LineWidth', 2);
            plot([i-dx, i+dx], [R_avg_avg_avg(i)+R_avg_avg_SEs(i), R_avg_avg_avg(i)+R_avg_avg_SEs(i)], '-k', 'LineWidth', 2);
        end;
        for j = 1:3
            plot([bx(j,1)-dx, bx(j,2)+dx], [by(j), by(j)], '-k', 'LineWidth', 2);
            plot([bx(j,1)-dx, bx(j,1)-dx], [by(j)+0.001, by(j)-dy], '-k', 'LineWidth', 2);
            plot([bx(j,2)+dx, bx(j,2)+dx], [by(j)+0.001, by(j)-dy], '-k', 'LineWidth', 2);
            text(mean(bx(j,:)), by(j), '*', 'FontSize', 28, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        end;
        axis([(1-1), (3+1), 0, 0.5]);
        set(gca,'Box','On');
        set(gca,'XTick', [1:3], 'XTickLabel', method);
        xlabel('decoding approach', 'FontSize', 16);
        ylabel('correlation coefficient', 'FontSize', 16);
        title('comparison', 'FontSize', 20); 
    end;
end;

cd('..');

% display colorbars
figure('Name', 'decoding accuracies', 'Color', [1 1 1], 'Position', [50 50 1200 1200]);
sp = [1, 4, 7, 2, 5, 8];

for i = 1:numel(sp)
    subplot(3,3,sp(i));
    if mod(sp(i),3) == 1, caxis([0, 2/3]); end;
    if mod(sp(i),3) == 2, caxis([0, 1/6]); end;
    axis off;
    colorbar('Location', 'SouthOutside');
end;