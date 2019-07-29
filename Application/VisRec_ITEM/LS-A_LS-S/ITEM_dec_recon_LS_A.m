function ITEM_dec_recon_LS_A(SPM, ROI, c, con, reg)
% _
% Decoding from Trials for Reconstruction (***least squares, all***)
% FORMAT ITEM_dec_recon(SPM, ROI, c, con, reg)
%     SPM - a structure specifying an estimated GLM
%     ROI - a filepath to a region of interest image
%     c   - a 1 x p vector with +1s indicating variables to reconstruct or
%           a q x p matrix with one +1 in each row indexing each variable
%     con - a string without spaces describing the variables (e.g. 'PE')
%     reg - a string without spaces describing the region (e.g. 'NAcc')
% 
% FORMAT ITEM_dec_recon(SPM, ROI, c, con, reg) performs decoding by an
% inverted transformation encoding model for reconstruction of selected
% variables indicated by c using the voxels included in ROI.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 30/11/2018, 09:10 (V0.1)
%  Last edit: 19/02/2019, 11:15 (Vn/a)


%=========================================================================%
% P R E P A R A T I O N                                                   %
%=========================================================================%

% Get SPM.mat if necessary
%-------------------------------------------------------------------------%
if nargin == 0
    SPM_mat = spm_select(1,'^SPM\.mat$','Select SPM.mat!');
    SPM_dir = fileparts(SPM_mat); load(SPM_mat);
    SPM.swd = SPM_dir;
    ITEM_dec_recon_LS_A(SPM);
    return
end;

% Set region of interest if necessary
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(ROI)
    ROI = strcat(SPM.swd,'/','mask.nii');
end;

% Set decoding contrast if necessary
%-------------------------------------------------------------------------%
if nargin < 3 || isempty(c)
    c = zeros(1,0);
    k = 0;
    while isempty(c) && k < numel(SPM.Sess(1).U(1))
    	k = k + 1;
        if ~strcmp(SPM.Sess(1).U(k).P(1).name,'none')
            c(1,1) = SPM.Sess(1).U(k).P(1).i(2);
        end;
    end;
    if isempty(c)
        c(1,1) = 1;
    end;
end;

% Set contrast name if necessary
%-------------------------------------------------------------------------%
if nargin < 4 || isempty(con)
    if size(c,1) > 1, d = sum(c,1); else d = c; end;
    con = '';
    ind = find(d);
    for l = 1:numel(ind)
        con = strcat(con,int2str(ind(l)));
        if l < numel(ind), con = strcat(con,','); end;
    end;
    clear d ind
end;

% Set region name if necessary
%-------------------------------------------------------------------------%
if nargin < 5 || isempty(reg)
    [fdir, reg, fext] = fileparts(ROI);
    clear fdir fext
end;

% Change to SPM.swd if specified
%-------------------------------------------------------------------------%
orig_dir = pwd;
try
    cd(SPM.swd);
catch
    SPM.swd = pwd;
end

% Get number of sessions
%-------------------------------------------------------------------------%
s = numel(SPM.Sess);

% Load GLM.mats in sub-directory
%-------------------------------------------------------------------------%
load(strcat(SPM.swd,'/','ITEM_est_1st_lvl_LS_A','/','GLM1.mat'));

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','ITEM_dec_recon_LS_A: load');
spm_progress_bar('Init',100,'Load trial-wise parameter estimates...','');

% Get ROI voxels
%-------------------------------------------------------------------------%
[roi_ind] = ITEM_load_ROI(GLM1, ROI);

% Load gamma estimates
%-------------------------------------------------------------------------%
t = sum(GLM1.tr);
v = numel(roi_ind);
d = ceil(t/100);
G = zeros(t,v);
for h = 1:numel(SPM.Sess)
    for k = 1:GLM1.tr(h)
        i = GLM1.Sess(h).t(k);
        g_str  = strcat(GLM1.swd,'/',GLM1.Vgamma(i).fname);
        g_hdr  = spm_vol(g_str);
        g_img  = spm_read_vols(g_hdr);
        G(i,:) = g_img(roi_ind);
        if mod(i,d) == 0, spm_progress_bar('Set',(i/t)*100); end;
    end;
end;
clear g_str g_hdr g_img

% Remove NaN voxels
%-------------------------------------------------------------------------%
v = v - sum(sum(isnan(G),1)>0);
G = G(:,sum(isnan(G),1)==0);

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');


%=========================================================================%
% E S T I M A T I O N   ( 1 ) :   I N V E R S E   M O D E L               %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','ITEM_dec_recon_LS_A: estimate (1)');

% Augment decoding contrast if necessary
%-------------------------------------------------------------------------%
if size(c,1) > 1, c = sum(c,1); end;
c = [c, zeros(1, GLM1.p(1)-numel(c))];
q = numel(find(c));

% Cycle through recording sessions
%-------------------------------------------------------------------------%
for h = 1:s
    
    % "data" - the T matrix
    %---------------------------------------------------------------------%
    Yh = GLM1.Sess(h).T(1:GLM1.t(h),1:GLM1.p(h));
    ITEM.Sess(h).Y = Yh;
    clear Yh
    
    % "design" - gamma estimates
    %---------------------------------------------------------------------%
    Xh = [G(GLM1.Sess(h).t(1:GLM1.t(h)),:), ones(GLM1.t(h),1)];
    ITEM.Sess(h).X = Xh;
    clear Xh
    
    % "covariance" - the U matrix
    %---------------------------------------------------------------------%
    Yh    = G(GLM1.Sess(h).t(1:GLM1.t(h)),:);
    Xh    = GLM1.Sess(h).T(1:GLM1.t(h),1:GLM1.p(h));
    Qh{1} = eye(GLM1.t(h));
    Qh{2} = GLM1.Sess(h).U(1:GLM1.t(h),1:GLM1.t(h));
    
    % Restricted maximum likelihood
    %---------------------------------------------------------------------%
    % [Vh, s2] = ITEM_GLM_ReML(Yh, Xh, Qh{1}, Qh{2}, sprintf('ITEM_dec_recon: ReML estimation for session %d',h));
    ITEM.Sess(h).V  = Qh{1}; % Qh{2}; % Vh;
    ITEM.Sess(h).s2 = [1 0]; % [0 1]; % s2;
    clear Yh Xh Qh Vh s2
    
end;


%=========================================================================%
% E S T I M A T I O N   ( 2 ) :   C R O S S - V A L I D A T I O N         %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','ITEM_dec_recon_LS_A: estimate (2)');

% Cycle through cross-validation folds
%-------------------------------------------------------------------------%
for g = 1:s
    
    % List sessions for this fold
    %---------------------------------------------------------------------%
    fold = 1:s;
    fold = fold(fold~=g);
    
    % Establish (in-sample) training data set
    %---------------------------------------------------------------------%
    Y_in = vertcat(ITEM.Sess(fold).Y);
    X_in = vertcat(ITEM.Sess(fold).X);
    V_in = blkdiag(ITEM.Sess(fold).V);
    
    % Establish (out-of-sample) test data set
    %---------------------------------------------------------------------%
    Y_out = ITEM.Sess(g).Y;
    X_out = ITEM.Sess(g).X;
    V_out = ITEM.Sess(g).V;
    W_out = sqrtm(inv(V_out));
    
    % Learn reconstruction weights from training data (concatenate sessions)
    %---------------------------------------------------------------------%
    [B_in, s2_in] = ITEM_GLM_MLE(Y_in, X_in, V_in, sprintf('Estimate weights that serve for reconstruction in session %d',g));
    
    % Apply reconstruction weights to test data
    %---------------------------------------------------------------------%
    Y_true  = Y_out(:,c==1);
    Y_pred  = W_out * X_out * B_in;
    Y_recon = Y_pred(:,c==1);
    
    % Calculate (out-of-sample) correlation coefficients
    %---------------------------------------------------------------------%
    CCg = diag(corr(Y_true,Y_recon))';
    ITEM.Sess(g).W  = W_out;
    ITEM.Sess(g).Yt = Y_true;
    ITEM.Sess(g).Yp = Y_pred;
    ITEM.Sess(g).Yr = Y_recon;
    ITEM.Sess(g).CC = CCg;
    clear CCg
    
end;


%=========================================================================%
% S A V E   R E S U L T S                                                 %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','ITEM_dec_recon_LS_A: save');

% Create target directory
%-------------------------------------------------------------------------%
ITEM.swd = strcat(SPM.swd,'/','ITEM_dec_recon_LS_A','/');
if ~exist(ITEM.swd,'dir'), mkdir(ITEM.swd); end;

% Complete ITEM structure
%-------------------------------------------------------------------------%
ITEM.ROI       = ROI;
ITEM.Recon.c   = c;
ITEM.Recon.con = con;
ITEM.Recon.reg = reg;

% Save ITEM structure
%-------------------------------------------------------------------------%
save(strcat(ITEM.swd,'ITEM_',con,'_',reg,'.mat'),'ITEM');

% Return to origin
%-------------------------------------------------------------------------%
cd(orig_dir);