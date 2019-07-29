% ITEM: Simulation
% _
% This script repeats the simulation study from Mumford et al. (2012) and
% extends it with an inverse transformed encoding model (ITEM) approach.
% 
% In the entire script, we use the following suffixes:
% - "*A": trial-wise design matrix X_S a.k.a.
%         Mumford's "least squares, all" (LS-A)
% - "*S": trial-based design matrix X_T a.k.a.
%         Mumford's "least squares, separate" (LS-S)
% - "*W": whitened trial-wise parameter estimates a.k.a.
%         Soch's "least squares, whitened" (LS-W)
% - "*T": (inverse) transformed encoding model a.k.a.
%         Soch's "least squares, transformed" (LS-T)
% 
% Author: Joram Soch
% E-Mail: joram.soch@bccn-berlin.de
% 
% Version History:
% - 29/11/2018, 12:30: finished simulation
% - 22/07/2019, 16:20: prepared for upload


clear
close all

%%% Step 1: set simulation parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n-> Step 1: ');

% set parameters
N   = 1000;                     % number of simulations
S   = 2;                        % number of sessions per test
t   = 60;                       % number of trials per session
bmu = [5, 3];                   % mean of the betas
bs2 = [0.5, 0.5];               % std  of the betas
es2 = [0.8; 1.6; 3];            % stds of the error (3)
t0  = 10;                       % stimulus onset
dt  = 2;                        % stimulus duration
TR  = 2;                        % repetition time
ISI = [0,4; 2,6; 4,8];          % inter-stimulus-intervals (3)
rho = 0.12;                     % temporal auto-correlation

% preallocate results
Sim    = struct([]);            % simulations structure
Res    = struct([]);            % sim results structure
TPR_bA = zeros(numel(es2),size(ISI,1));
TPR_bS = zeros(numel(es2),size(ISI,1));
TPR_bW = zeros(numel(es2),size(ISI,1));
TPR_bT = zeros(numel(es2),size(ISI,1));

fprintf('end.');


%%% Step 2: generate designs and data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n-> Step 2: ');

% for each noise level
for g = 1:numel(es2)
    
    fprintf('%1.1f: ', es2(g));
    
    % for each ISI range
    for h = 1:size(ISI,1)
        
        fprintf('[%d,%d], ', ISI(h,1), ISI(h,2));
        
        %%% Step 2a: set design fundamentals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % calculate number of scans needed for ISI
        n = ceil((t0 + t*dt + (t-1)*((ISI(h,1)+ISI(h,2))/2) + 5*t0)/TR);
        z = zeros(n,1);
        Sim(g,h).base.n = n;
        
        % create temporal auto-correlation accordingly
        V = toeplitz(rho.^[0:1:(n-1)]);
        W = sqrtm(inv(V));
        Sim(g,h).base.V = V;
        Sim(g,h).base.W = W;
        
        % multiply with variance to get covariance matrix
        s2V = es2(g)^2 * V;
        Sim(g,h).base.s2V = s2V;
        
        % configure settings for HRF convolution
        settings.n  = n;
        settings.TR = TR;
        
        %%% Step 2b: create design matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % for each session
        for j = 1:S
        
            % sample trial types
            ttm = [[1*ones(t/2,1); 2*ones(t/2,1)], rand(t,1)];
            ttm = sortrows(ttm, 2);
            tt  = ttm(:,1);
            Sim(g,h).des(j).tt = tt;

            % sample trial onsets
            isi = MD_unirnd(ISI(h,1), ISI(h,2), t-1)';
            ons = [t0, t0 + cumsum(isi) + dt*[1:(t-1)]]';
            dur = dt*ones(t,1);
            Sim(g,h).des(j).isi = isi';
            Sim(g,h).des(j).ons = ons;
            Sim(g,h).des(j).dur = dur;

            % create design matrix X_S (Mumford: LS-A)
            clear names onsets durations
            for k = 1:t
                names{k}     = strcat('trl-',MF_int2str0(k,2));
                onsets{k}    = ons(k);
                durations{k} = dur(k);
            end;
            [X, L] = ITEM_get_des_mat(names, onsets, durations, [], [], [], settings);
            Sim(g,h).des(j).DM.XS = X;

            % create design matrices X_T (Mumford: LS-S)
            X_t = X;
            for k = 1:t
                T = zeros(t,2);
                T(       k,1) = 1;
                T([1:t]~=k,2) = 1;
                X = X_t * T;
                Sim(g,h).des(j).DM.XT{k} = X;
            end;

            % create design matrix T (Soch: LS-W/LS-T)
            X = Sim(g,h).des(j).DM.XS;
            U = ((W*X)'*(W*X))^(-1);
            T = zeros(t,2);
            for k = 1:2
                T(tt==k,k) = 1;
            end;
            X = X*T;
            Sim(g,h).des(j).DM.X = X;
            Sim(g,h).des(j).DM.T = T;
            Sim(g,h).des(j).DM.U = U;
        
        end;
        
        %%% Step 2c: sample data matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % for each session
        for j = 1:S
            
            % sample parameters
            tt = Sim(g,h).des(j).tt;
            bA = MD_matnrnd(bmu(1)*ones(t/2, N), bs2(1)*eye(t/2), eye(N));
            bB = MD_matnrnd(bmu(2)*ones(t/2, N), bs2(2)*eye(t/2), eye(N));
          % bA = normrnd(bmu(1),bs2(1),[(t/2) N]);
          % bB = normrnd(bmu(2),bs2(2),[(t/2) N]);
            b  = zeros(t,N);
            b(tt==1,:) = bA;
            b(tt==2,:) = bB;
            Sim(g,h).data(j).b = b;
            clear bA bB

            % sample data
            X = Sim(g,h).des(j).DM.XS;
            e = MD_matnrnd(zeros(n,N), s2V, eye(N));
          % e = mvnrnd(z, s2V, N)';
            y = X*b + e;
            Sim(g,h).data(j).y = y;
            clear y e

        end;
        
    end;
    
    fprintf('done; ');
    
end;

fprintf('end.');


%%% Step 3: estimate models and test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n-> Step 3: ');

% for each noise level
for g = 1:numel(es2)
    
    fprintf('%1.1f: ', es2(g));
    
    % for each ISI range
    for h = 1:size(ISI,1)
        
        fprintf('[%d,%d]: \n', ISI(h,1), ISI(h,2));
        
        %%% Step 3a: estimate model parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % for each simulation
        for i = 1:N
            
            % for each session
            for j = 1:S
                
                % filter data and design
                % not necessary as not part of the ground truth
                
                % get data and whitening matrix
                y = Sim(g,h).data(j).y(:,i);
                W = Sim(g,h).base.W;
                
                % estimate using trial-wise design matrix
                Wy = W*y;
                WX = W*Sim(g,h).des(j).DM.XS;
                bA = (WX'*WX)^(-1) * WX'*Wy;
                Sim(g,h).est(j).bA(:,i) = bA;
                
                % estimate using trial-based design matrices
                bS = zeros(t,1);
                for k = 1:t
                    WX = W*Sim(g,h).des(j).DM.XT{k};
                    bk = (WX'*WX)^(-1) * WX'*Wy;
                    bS(k) = bk(1);
                end;
                Sim(g,h).est(j).bS(:,i) = bS;
                
                % whitening of trial-wise parameter estimates
                U  = Sim(g,h).des(j).DM.U;
                WU = sqrtm(inv(U));
                bW = WU*bA;
                Sim(g,h).est(j).bW(:,i) = bW;
                
            end;
            
        end;
        
        %%% Step 3b: restricted maximum likelihood %%%%%%%%%%%%%%%%%%%%%%%%
        
        % for each session
        for j = 1:S

            % prepare ReML analysis
            BA = Sim(g,h).est(j).bA;
            YY = (1/N) * (BA*BA');
            X  = Sim(g,h).des(j).DM.T;
            Q{1} = eye(t);
            Q{2} = Sim(g,h).des(j).DM.U;

            % perform ReML analysis
            [V, s2] = spm_reml(YY, X, Q);
            Sim(g,h).des(j).DM.Sg = V;
            Sim(g,h).des(j).DM.s2 = full(s2)';
            
        end;
        
        % for each simulation
        for i = 1:N
            
            % for each session
            for j = 1:S
                
                %%% Step 3c: perform statistical tests %%%%%%%%%%%%%%%%%%%%
                
                % two-sample t-test on trial-wise parameter estimates
                tt = Sim(g,h).des(j).tt;
                bA = Sim(g,h).est(j).bA(:,i);
                [H, p] = ttest2(bA(tt==1), bA(tt==2));
                Sim(g,h).test.hA(j,i) = H;
                Sim(g,h).test.pA(j,i) = p;
                
                % two-sample t-test on trial-based parameter values
                bS = Sim(g,h).est(j).bS(:,i);
                [H, p] = ttest2(bS(tt==1), bS(tt==2));
                Sim(g,h).test.hS(j,i) = H;
                Sim(g,h).test.pS(j,i) = p;
                
                % two-sample t-test on whitened parameter estimates
                bW = Sim(g,h).est(j).bW(:,i);
                [H, p] = ttest2(bW(tt==1), bW(tt==2));
                Sim(g,h).test.hW(j,i) = H;
                Sim(g,h).test.pW(j,i) = p;
                
                % two-sample t-test using transformation encoding model
                y = Sim(g,h).est(j).bA(:,i);
                X = Sim(g,h).des(j).DM.T;
                V = Sim(g,h).des(j).DM.Sg;
                [H, p] = ITEM_GLM_con(y, X, V, [1 -1]', 'F', 0.05);
                Sim(g,h).test.hT(j,i) = H;
                Sim(g,h).test.pT(j,i) = p;
                
                %%% Step 3d: decode/classify trials %%%%%%%%%%%%%%%%%%%%%%%
                
                % logistic regression using trial-wise parameter estimates
                y_train = Sim(g,h).des([1:S]~=j).tt;
                y_test  = Sim(g,h).des(j).tt;
                X_train = Sim(g,h).est([1:S]~=j).bA(:,i);
                X_test  = Sim(g,h).est(j).bA(:,i);
                b_train = mnrfit(X_train,y_train);
                o_pred  = exp(X_test*b_train(2) + b_train(1));
                y_pred  = zeros(size(y_test));
                y_pred(o_pred>1) = min(y_train);
                y_pred(o_pred<1) = max(y_train);
                a = sum(y_test==y_pred)/t;
                Sim(g,h).pred.aA(j,i) = a;
                
                % logistic regression using trial-based parameter values
                X_train = Sim(g,h).est([1:S]~=j).bS(:,i);
                X_test  = Sim(g,h).est(j).bS(:,i);
                b_train = mnrfit(X_train,y_train);
                o_pred  = exp(X_test*b_train(2) + b_train(1));
                y_pred(o_pred>1) = min(y_train);
                y_pred(o_pred<1) = max(y_train);
                a = sum(y_test==y_pred)/t;
                Sim(g,h).pred.aS(j,i) = a;
                
                % logistic regression using whitened parameter estimates
                X_train = Sim(g,h).est([1:S]~=j).bW(:,i);
                X_test  = Sim(g,h).est(j).bW(:,i);
                b_train = mnrfit(X_train,y_train);
                o_pred  = exp(X_test*b_train(2) + b_train(1));
                y_pred(o_pred>1) = min(y_train);
                y_pred(o_pred<1) = max(y_train);
                a = sum(y_test==y_pred)/t;
                Sim(g,h).pred.aW(j,i) = a;
                
                % decoding using inverted transformation encoding model
                Y_train = Sim(g,h).des([1:S]~=j).DM.T;
                Y_test  = Sim(g,h).des(j).DM.T;
                X_train = [Sim(g,h).est([1:S]~=j).bA(:,i), ones(t,1)];
                X_test  = [Sim(g,h).est(j).bA(:,i), ones(t,1)];
                V_train = Sim(g,h).des([1:S]~=j).DM.Sg;
                V_test  = Sim(g,h).des(j).DM.Sg;
                P_train = inv(V_train);
                W_test  = sqrtm(inv(V_test));
                b_train = (X_train'*P_train*X_train)^(-1) * X_train'*P_train*Y_train;
                Y_rec   = W_test * X_test * b_train;
                Y_pred  = [(Y_rec(:,1)>Y_rec(:,2)), (Y_rec(:,2)>Y_rec(:,1))];
                a = sum(Y_test(:,1)==Y_pred(:,1))/t;
                Sim(g,h).pred.aT(j,i) = a;
                
            end;
            
        end;
        
    end;
    
    fprintf('done; ');
    
end;

fprintf('end.');


%%% Step 4: compute results summaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n-> Step 4: ');

% for each noise level
for g = 1:numel(es2)
    
    % for each ISI range
    for h = 1:size(ISI,1)
        
        % for each simulation
        for i = 1:N
        
            % for each session
            for j = 1:S

                % trial types
                tt = Sim(g,h).des(j).tt;
                
                % mean squared errors
                % for k = 1:2
                %     Res(g,h).MSE(j).eA(k,i) = mean((Sim(g,h).est(j).bA(tt==k,i) - Sim(g,h).data(j).b(tt==k,i)).^2);
                %     Res(g,h).MSE(j).eS(k,i) = mean((Sim(g,h).est(j).bS(tt==k,i) - Sim(g,h).data(j).b(tt==k,i)).^2);
                %     Res(g,h).MSE(j).eW(k,i) = mean((Sim(g,h).est(j).bW(tt==k,i) - Sim(g,h).data(j).b(tt==k,i)).^2);
                % end;

                % mean squared errors
                Res(g,h).MSE(j).eA(1,i) = mean((Sim(g,h).est(j).bA(:,i) - Sim(g,h).data(j).b(:,i)).^2);
                Res(g,h).MSE(j).eS(1,i) = mean((Sim(g,h).est(j).bS(:,i) - Sim(g,h).data(j).b(:,i)).^2);
                Res(g,h).MSE(j).eW(1,i) = mean((Sim(g,h).est(j).bW(:,i) - Sim(g,h).data(j).b(:,i)).^2);

                % correlations
                for k = 1:2
                    Res(g,h).corr(j).rA(k,i) = corr(Sim(g,h).data(j).b(tt==k,i), Sim(g,h).est(j).bA(tt==k,i));
                    Res(g,h).corr(j).rS(k,i) = corr(Sim(g,h).data(j).b(tt==k,i), Sim(g,h).est(j).bS(tt==k,i));
                    Res(g,h).corr(j).rW(k,i) = corr(Sim(g,h).data(j).b(tt==k,i), Sim(g,h).est(j).bW(tt==k,i));
                end;
                
                % variances
                for k = 1:2
                    Res(g,h).var(j).s2A(k,i) = var(Sim(g,h).est(j).bA(tt==k,i));
                    Res(g,h).var(j).s2S(k,i) = var(Sim(g,h).est(j).bS(tt==k,i));
                    Res(g,h).var(j).s2W(k,i) = var(Sim(g,h).est(j).bW(tt==k,i));
                end;
                
                % auto-correlations
                for k = 1:3
                    Res(g,h).acorr(j).rA(k,i) = corr(Sim(g,h).est(j).bA(1:end-k,i), Sim(g,h).est(j).bA(1+k:end,i));
                    Res(g,h).acorr(j).rS(k,i) = corr(Sim(g,h).est(j).bS(1:end-k,i), Sim(g,h).est(j).bS(1+k:end,i));
                    Res(g,h).acorr(j).rW(k,i) = corr(Sim(g,h).est(j).bW(1:end-k,i), Sim(g,h).est(j).bW(1+k:end,i));
                end;
                
            end;
        
        end;
        
        % average across sessions
        Res(g,h).all.MSE  = [mean(vertcat(Res(g,h).MSE.eA));  mean(vertcat(Res(g,h).MSE.eS));  mean(vertcat(Res(g,h).MSE.eW))];
        Res(g,h).all.corr = [mean(vertcat(Res(g,h).corr.rA)); mean(vertcat(Res(g,h).corr.rS)); mean(vertcat(Res(g,h).corr.rW))];
        Res(g,h).all.var  = [mean(vertcat(Res(g,h).var.s2A)); mean(vertcat(Res(g,h).var.s2S)); mean(vertcat(Res(g,h).var.s2W))];
        for k = 1:3
            Res(g,h).all.acorr{k} = [mean([Res(g,h).acorr(1).rA(k,:); Res(g,h).acorr(2).rA(k,:)]);
                                     mean([Res(g,h).acorr(1).rS(k,:); Res(g,h).acorr(2).rS(k,:)]);
                                     mean([Res(g,h).acorr(1).rW(k,:); Res(g,h).acorr(2).rW(k,:)])];
        end;
        
        % test performances
        TPR_bA(g,h) = sum(sum(Sim(g,h).test.hA))/(S*N);
        TPR_bS(g,h) = sum(sum(Sim(g,h).test.hS))/(S*N);
        TPR_bW(g,h) = sum(sum(Sim(g,h).test.hW))/(S*N);
        TPR_bT(g,h) = sum(sum(Sim(g,h).test.hT))/(S*N);
        
        % decoding accuracies
        Res(g,h).DA = [mean(Sim(g,h).pred.aA); mean(Sim(g,h).pred.aS); mean(Sim(g,h).pred.aW); mean(Sim(g,h).pred.aT)];
        
    end;
    
end;

save('Simulation.mat', 'Sim', 'Res', 'TPR_bA', 'TPR_bS', 'TPR_bW', 'TPR_bT');

fprintf('end.');


%%% Step 5: plot simulation results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n');

% plot design matrices
figure('Name', 'design matrices', 'Color', [1 1 1], 'Position', [50 50 900 900]);

subplot(5,3,[1 4 7]);
imagesc(Sim(1,1).des(1).DM.X); axis off;
title('X_{ } [n x p]', 'FontSize', 16);

subplot(5,3,[2 5 8]);
imagesc(Sim(1,1).des(1).DM.XS); axis off;
title('X_t [n x t]', 'FontSize', 16);

subplot(5,3,[3 6 9]);
imagesc(Sim(1,1).des(1).DM.T); axis off;
title('T_{ } [t x p]', 'FontSize', 16);

subplot(5,3,[11 14]);
imagesc(Sim(1,1).des(1).DM.U); axis off; axis square;
title('U = (X_{t}^{T}V^{-1}X_{t})^{-1}', 'FontSize', 16);

% plot linear models
figure('Name', 'linear models', 'Color', [1 1 1], 'Position', [50 50 900 900]);

% standard model
subplot(3,4,1);
axis off
text(1/2, 1/2, 'standard model', 'FontSize', 16, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');

subplot(3,4,2);
imagesc(Sim(1,1).data(1).y(:,1)); axis off;
title('y', 'FontSize', 16);

subplot(3,4,3);
imagesc(Sim(1,1).des(1).DM.X); axis off;
title('X', 'FontSize', 16);

subplot(3,4,4);
imagesc(Sim(1,1).base.V); axis off; axis square;
title('V', 'FontSize', 16);

% first-level model
subplot(3,4,5);
axis off
text(1/2, 1/2, '"first-level" model', 'FontSize', 16, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');

subplot(3,4,6);
imagesc(Sim(1,1).data(1).y(:,1)); axis off;
title('y', 'FontSize', 16);

subplot(3,4,7);
imagesc(Sim(1,1).des(1).DM.XS); axis off;
title('X_t', 'FontSize', 16);

subplot(3,4,8);
imagesc(Sim(1,1).base.V); axis off; axis square;
title('V', 'FontSize', 16);

% second-level model
subplot(3,4,9);
axis off
text(1/2, 1/2, '"second-level" model', 'FontSize', 16, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');

subplot(3,4,10);
imagesc(Sim(1,1).est(1).bA(:,1)); axis off;
title('\gamma', 'FontSize', 16);

subplot(3,4,11);
imagesc(Sim(1,1).des(1).DM.T); axis off;
title('T', 'FontSize', 16);

subplot(3,4,12);
imagesc(Sim(1,1).des(1).DM.U); axis off; axis square;
title('U', 'FontSize', 16);

% plot mean squared errors
figure('Name', 'mean squared errors', 'Color', [1 1 1], 'Position', [50 50 1000 1000]);

for g = 1:numel(es2)
    for h = 1:size(ISI,1)
        subplot(numel(es2), size(ISI,1), (g-1)*size(ISI,1)+h);
        hBP = boxplot(Res(g,h).all.MSE','Positions',[1:3],'Width',2/3,'Colors','rbc','Symbol','+k','Labels',{'LS-A' 'LS-S' 'LS-W'});
        xlim([(1-0.5), (3+0.5)]);
        set(hBP,'LineWidth',2);
        set(gca,'Box','On');
        if g == numel(es2)
            xlabel(['t_{isi} \sim', sprintf(' U(%d,%d)', ISI(h,1), ISI(h,2))], 'FontSize', 16);
        end;
        if h == 1
            ylabel(['\sigma^2 =', sprintf(' %1.1f', es2(g))], 'FontSize', 16);
        end;
        if g == 1 & h == floor(mean([1 size(ISI,1)]))
            title('Mean Squared Errors', 'FontSize', 20);
        end;
    end;
end;

% plot correlations
figure('Name', 'correlations', 'Color', [1 1 1], 'Position', [50 50 1000 1000]);

for g = 1:numel(es2)
    for h = 1:size(ISI,1)
        subplot(numel(es2), size(ISI,1), (g-1)*size(ISI,1)+h);
        hBP = boxplot(Res(g,h).all.corr','Positions',[1:3],'Width',2/3,'Colors','rbc','Symbol','+k','Labels',{'LS-A' 'LS-S' 'LS-W'});
        xlim([(1-0.5), (3+0.5)]);
        ylim([0 1]);
        set(hBP,'LineWidth',2);
        set(gca,'Box','On');
        if g == numel(es2)
            xlabel(['t_{isi} \sim', sprintf(' U(%d,%d)', ISI(h,1), ISI(h,2))], 'FontSize', 16);
        end;
        if h == 1
            ylabel(['\sigma^2 =', sprintf(' %1.1f', es2(g))], 'FontSize', 16);
        end;
        if g == 1 & h == floor(mean([1 size(ISI,1)]))
            title('Correlations', 'FontSize', 20);
        end;
    end;
end;

% plot variances
figure('Name', 'variances', 'Color', [1 1 1], 'Position', [50 50 1000 1000]);

for g = 1:numel(es2)
    for h = 1:size(ISI,1)
        subplot(numel(es2), size(ISI,1), (g-1)*size(ISI,1)+h);
        hBP = boxplot(Res(g,h).all.var','Positions',[1:3],'Width',2/3,'Colors','rbc','Symbol','+k','Labels',{'LS-A' 'LS-S' 'LS-W'});
        xlim([(1-0.5), (3+0.5)]);
        set(hBP,'LineWidth',2);
        set(gca,'Box','On');
        if g == numel(es2)
            xlabel(['t_{isi} \sim', sprintf(' U(%d,%d)', ISI(h,1), ISI(h,2))], 'FontSize', 16);
        end;
        if h == 1
            ylabel(['\sigma^2 =', sprintf(' %1.1f', es2(g))], 'FontSize', 16);
        end;
        if g == 1 & h == floor(mean([1 size(ISI,1)]))
            title('Variances', 'FontSize', 20);
        end;
    end;
end;

% plot auto-correlations
figure('Name', 'auto-correlations', 'Color', [1 1 1], 'Position', [50 50 1000 1000]);

for g = 1:numel(es2)
    for h = 1:size(ISI,1)
        sp = subplot(numel(es2), size(ISI,1), (g-1)*size(ISI,1)+h);
        hold on;
        for k = 1:3
            hBP = boxplot(Res(g,h).all.acorr{k}','Positions',k+[-1/4, 0, 1/4],'Width',1/6,'Colors','rbc','Symbol','+k','Labels',{'' '' ''});
            set(hBP,'LineWidth',2);
            text(k, -1.01, num2str(k), 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Top');
        end;
        xlim([(1-0.5), (3+0.5)]);
        ylim([-1 +1]);
        set(sp,'Box','On');
        if g == numel(es2)
            xlabel(['t_{isi} \sim', sprintf(' U(%d,%d)', ISI(h,1), ISI(h,2))], 'FontSize', 16);
        end;
        if h == 1
            ylabel(['\sigma^2 =', sprintf(' %1.1f', es2(g))], 'FontSize', 16);
        end;
        if g == 1 & h == floor(mean([1 size(ISI,1)]))
            title('Auto-Correlations', 'FontSize', 20);
            xlabel('lag', 'FontSize', 20);
        end;
    end;
end;


% plot test performances
figure('Name', 'test performances', 'Color', [1 1 1], 'Position', [50 50 1000 1000]);

for g = 1:numel(es2)
    for h = 1:size(ISI,1)
        subplot(numel(es2), size(ISI,1), (g-1)*size(ISI,1)+h);
        hold on;
        bar(1, TPR_bA(g,h), 'r');
        bar(2, TPR_bS(g,h), 'b');
        bar(3, TPR_bW(g,h), 'c');
        bar(4, TPR_bT(g,h), 'g');
        if diff(bmu) ~= 0
            axis([(1-0.5), (4+0.5), 0.5, 1.01]);
        else
            axis([(1-0.5), (4+0.5), -0.01, 0.5]);
        end;
        set(gca,'Box','On');
        set(gca,'XTick',[1:4],'XTickLabel',{'LS-A' 'LS-S' 'LS-W' 'LS-T'});
        if g == numel(es2)
            xlabel(['t_{isi} \sim', sprintf(' U(%d,%d)', ISI(h,1), ISI(h,2))], 'FontSize', 16);
        end;
        if h == 1
            ylabel(['\sigma^2 =', sprintf(' %1.1f', es2(g))], 'FontSize', 16);
        end;
        if g == 1 & h == floor(mean([1 size(ISI,1)]))
            if diff(bmu) ~= 0
                title('True Positive Rates', 'FontSize', 20);
            else
                title('False Positive Rates', 'FontSize', 20);
            end;
        end;
    end;
end;

% plot decoding accuracies
figure('Name', 'decoding accuracies', 'Color', [1 1 1], 'Position', [50 50 1000 1000]);

for g = 1:numel(es2)
    for h = 1:size(ISI,1)
        subplot(numel(es2), size(ISI,1), (g-1)*size(ISI,1)+h);
        hBP = boxplot(Res(g,h).DA','Positions',[1:4],'Width',2/3,'Colors','rbcg','Symbol','+k','Labels',{'LS-A' 'LS-S' 'LS-W' 'LS-T'});
        set(hBP,'LineWidth',2);
        if diff(bmu) ~= 0
            axis([(1-0.5), (4+0.5), 0.5, 1]);
        else
            axis([(1-0.5), (4+0.5), 0.25, 0.75]);
        end;
        set(gca,'Box','On');
        if g == numel(es2)
            xlabel(['t_{isi} \sim', sprintf(' U(%d,%d)', ISI(h,1), ISI(h,2))], 'FontSize', 16);
        end;
        if h == 1
            ylabel(['\sigma^2 =', sprintf(' %1.1f', es2(g))], 'FontSize', 16);
        end;
        if g == 1 & h == floor(mean([1 size(ISI,1)]))
            title('Decoding Accuracies', 'FontSize', 20);
        end;
    end;
end;