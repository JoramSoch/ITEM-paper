% ITEM Paper, Figure 5 (old Figure 4)

clear
close all

% load decoding accuracies
load Simulation.mat
bmu = [5, 3];                   % mean of the betas
bs2 = [0.5, 0.5];               % std  of the betas
es2 = [0.8; 1.6; 3];            % stds of the error (3)
ISI = [0,4; 2,6; 4,8];          % inter-stimulus-intervals (3)

% display decoding accuracies
figure('Name', 'decoding accuracies', 'Color', [1 1 1], 'Position', [50 50 1000 1000]);

for g = 1:numel(es2)
    for h = 1:size(ISI,1)
        subplot(numel(es2), size(ISI,1), (g-1)*size(ISI,1)+h);
        DAs = Res(g,h).DA([1, 2, 4],:);
        hBP = boxplot(DAs','Positions',[1:3],'Width',2/3,'Colors','rbg','Symbol','+k','Labels',{'LS-A' 'LS-S' 'LS-T'});
        axis([(1-0.5), (3+0.5), 0.5, 1]);
        set(hBP,'LineWidth',2);
        set(gca,'Box','On');
        if g == numel(es2)
            xlabel(['t_{isi} \sim', sprintf(' U(%d,%d)', ISI(h,1), ISI(h,2))], 'FontSize', 16);
        end;
        if h == 1
            ylabel(['\sigma^2 =', sprintf(' %1.1f', es2(g))], 'FontSize', 16);
        end;
    end;
end;