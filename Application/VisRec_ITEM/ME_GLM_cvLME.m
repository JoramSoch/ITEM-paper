function [cvLME, oosLME] = ME_GLM_cvLME(GLM)
% _
% Cross-Validated Log Model Evidence for General Linear Model
% FORMAT [cvLME, oosLME] = ME_GLM_cvLME(GLM)
% 
%     GLM    - a  1 x s structure with the following fields:
%     o Y    - an n x v data matrix of v time series with n data points
%     o X    - an n x p design matrix of p regressors with n data points
%     o P    - an n x n precision matrix embodying covariance assumptions
% 
%     cvLME  - a  1 x v vector of cross-validated log model evidences
%     oosLME - an s x v matrix of out-of-sample log model evidences
% 
% FORMAT [cvLME, oosLME] = ME_GLM_cvLME(GLM) calculates the cross-validated
% log model evidence for a GLM describing s subsets with data matrix Y,
% design matrix X, precision matrix P and normal-gamma distributed priors
% for regression coefficients and residual variance.
% 
% Further information:
%     help ME_GLM_NG
%     help ME_GLM_NG_LME
%     help ME_GLM_NG_AnC
% 
% References:
% [1] Soch J, Haynes JD, Allefeld C (2016): "How to avoid mismodelling in
%     GLM-based fMRI data analysis: cross-validated Bayesian model selection".
%     NeuroImage, vol. 141, pp. 469–489.
% [2] Soch J, Meyer AP, Haynes JD, Allefeld C (2017): "How to improve parameter estimates in 
%     GLM-based fMRI data analysis: cross-validated Bayesian model averaging".
%     NeuroImage, vol. 158, pp. 186-195.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 21/02/2018, 14:20 (V1.2/V18)
%  Last edit: 21/02/2018, 14:20 (V1.2/V18)


% Get model dimensions
%-------------------------------------------------------------------------%
s = numel(GLM);                 % number of data subsets
v = size(GLM(1).Y,2);           % number of time series
n = size(GLM(1).Y,1);           % number of data points
p = size(GLM(1).X,2);           % number of parameters

% Set non-informative priors
%-------------------------------------------------------------------------%
m0_ni = zeros(p,1);             % flat Gaussian
L0_ni = exp(-23)*eye(p);
a0_ni = 0;                      % Jeffrey's prior
b0_ni = 0;

% Estimate posteriors from all data
%-------------------------------------------------------------------------%
Y = vertcat(GLM(1:s).Y);
X = vertcat(GLM(1:s).X);
P = blkdiag(GLM(1:s).P);
[mn_all, Ln_all, an_all, bn_all] = ME_GLM_NG(Y, X, P, m0_ni, L0_ni, a0_ni, b0_ni, sprintf('Estimate posteriors over all sessions 1-%d',s));
clear Y X P

% Preallocate out-of-sample LMEs
%-------------------------------------------------------------------------%
oosLME = zeros(s,v);

% Cycle through cross-validation folds
%-------------------------------------------------------------------------%
for i = 1:s
    
    % List sessions for this fold
    %---------------------------------------------------------------------%
    fold = 1:s;
    fold = fold(fold~=i);

    % Estimate posteriors for this fold
    %---------------------------------------------------------------------%
    Yf = vertcat(GLM(fold).Y);
    Xf = vertcat(GLM(fold).X);
    Pf = blkdiag(GLM(fold).P);
    [mn, Ln, an, bn] = ME_GLM_NG(Yf, Xf, Pf, m0_ni, L0_ni, a0_ni, b0_ni, sprintf('Estimate posteriors that serve as priors for session %d',i));
    clear Yf Xf Pf
    
    % Calculate out-of-sample LMEs
    %---------------------------------------------------------------------%
    oosLME(i,:) = ME_GLM_NG_LME(GLM(i).P, Ln, an, bn, Ln_all, an_all, bn_all);
    
end;

% Calculate cross-validated LMEs
%-------------------------------------------------------------------------%
cvLME = sum(oosLME,1);