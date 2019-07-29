function MLG = ME_MGLM_oosinv(GLM, train, test)
% _
% Out-of-Sample Model Inversion for Multivariate General Linear Model
% FORMAT MLG = ME_MGLM_oosinv(GLM, train, test)
% 
%     GLM   - a  1 x s structure with the following fields:
%     o Y   - an n x v data matrix of v time series with n data points
%     o X   - an n x p design matrix of p regressors with n data points
%     o V   - an n x n covariance matrix embodying covariance assumptions
%     o W   - an n x n whitening matrix for removing serial dependencies
%     train - a  1 x (s-1) vector indicating training subsets
%     test  - an integer <= s indicating the test subset
% 
%     MLG   - a  1 x s structure with the following fields:
%     o B   - a  p x v parameter matrix, the estimated activation pattern
%     o Sx  - a  p x p covariance matrix, the empirical design covariance
%     o Sy  - a  v x v covariance matrix, the empirical data covariance
%     o W   - a  v x p weight matrix, the estimated extraction filter
%     o X   - an n x p design matrix, the reconstructed design matrix
% 
% FORMAT MLG = ME_MGLM_oosinv(GLM, train, test) estimates session-wise
% parameters of a multivariate GLM describing s subsets with data matrix Y,
% design matrix X, covariance matrix V and performs an out-of-sample model
% inversion to reconstruct the design matrix of the left-out session.
% 
% References:
% [1] Haufe S, Meinecke F, Görgen K et al. (2014): "On the interpretation
%     of weight vectors of linear models in multivariate neuroimaging".
%     NeuroImage, vol. 87, pp. 96-110.
% [2] Soch J, Haynes JD, Allefeld C (2016): "How to avoid mismodelling in
%     GLM-based fMRI data analysis: cross-validated Bayesian model selection".
%     NeuroImage, vol. 141, pp. 469–489.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 21/02/2018, 15:25 (V1.2/V18)
%  Last edit: 11/04/2018, 13:30 (V1.2/V18)


% Get model dimensions
%-------------------------------------------------------------------------%
s = numel(GLM);                 % number of data subsets
v = size(GLM(1).Y,2);           % number of time series
n = size(GLM(1).Y,1);           % number of data points
p = size(GLM(1).X,2);           % number of parameters

% Preallocate inverted GLM
%-------------------------------------------------------------------------%
MLG = struct([]);

% Estimate session-wise parameters
%-------------------------------------------------------------------------%
for j = train
    
    % Get data, design, covariance, precision
    %---------------------------------------------------------------------%
    Y = GLM(j).Y;
    X = GLM(j).X;
    V = GLM(j).V;
    P = inv(V);
    
    % Get beta estimates, empirical covariances
    %---------------------------------------------------------------------%
    MLG(j).B  = (X'*P*X)^-1 * (X'*P*Y);
    MLG(j).Sx = 1/(n-1) * (X'*P*X);
    MLG(j).Sy = 1/(n-1) * (Y'*P*Y);
    
end;

% Average empirical covariances
%---------------------------------------------------------------------%
A  = zeros(v,p);
Sx = zeros(p,p);
Sy = zeros(v,v);
for j = train
    A  = A  + MLG(j).B';
    Sx = Sx + MLG(j).Sx;
    Sy = Sy + MLG(j).Sy;
end;
A  = 1/(s-1) * A;
Sx = 1/(s-1) * Sx;
Sy = 1/(s-1) * Sy;
i  = test;

% Transform into weight matrix
%---------------------------------------------------------------------%
MLG(i).W = pinv(Sy) * A * Sx;
GLM(i).W = sqrtm(inv(GLM(i).V));

% Reconstruct design matrix
%---------------------------------------------------------------------%
MLG(i).X = GLM(i).W * GLM(i).Y * MLG(i).W;