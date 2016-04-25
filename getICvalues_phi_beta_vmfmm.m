function valIC = getICvalues_phi_beta_vmfmm(X, param, labels, cp)
% Construct a hierarchy of vMF mixture models by exploiting Bregman
% Divergence and compute the scores of different information criteria
% See Section 4.3 of [1] for details

% INPUT:
% X: feature vectors (N x 3)
% param: parameters of a vMF mixture model
% labels: assigned labels to the observations w.r.t. the mixture model
% parameters
% cp: conditional probability

% OUTPUT
% valIC: scores of different information criteria

% Reference:
% [1] Hasnat et al., Model-based hierarchical clustering with Bregman 
% divergences and Fishers mixture model: application to depth image analysis. 
% Statistics and Computing, 1-20, 2015.

% Author: Md Abul HASNAT

N = size(X,1);

% compte beta
beta_aic = (log(3) - log(log(log(N))))/log(N);
beta_bic = (log(log(N)) - log(log(log(N))))/log(N);
beta_min = (log(log(N)) / (log(N)));
beta_max = 1 - beta_min;
beta_range = [0 beta_aic beta_bic beta_min:(beta_max-beta_min)/20:beta_max];

% compute C_phi_beta(k)
C_phi_beta_k_aic      = (N^beta_aic) * log(log(N));
C_phi_beta_k_bic      = (N^beta_bic) * log(log(N));
C_phi_beta_k_beta_min = (N^beta_min) * log(log(N));
C_phi_beta_k_beta_max = (N^beta_max) * log(log(N));
C_phi_beta_k_beta_range = (N.^beta_range) .* log(log(N));

% Get the information
mu = param.mu;
kappa = param.kappa;
alpha = param.weight;

k = length(alpha);
p = size(mu,2);

% compute P(k)
P_k = (k*(p + 1)) + k - 1;

%
%% Compute log likelihood for the optimal label
% Compute new log likelihood
logWeight   = log(alpha);
logNormTerm = log(kappa) - log(4*pi*sinh(kappa));
logExpTerm  = bsxfun(@times, kappa,  (mu * X')');

logClassCondLiklihood = bsxfun(@plus, logWeight + logNormTerm , logExpTerm);

% log likelihood value
ClassCondLiklihood    = exp(logClassCondLiklihood);
llh = sum(log(sum(ClassCondLiklihood,2)));

cp = cp + 0.00001;
lcp = log(cp);

sumLCPOp = 0;     % Sum of log conditional probability
for tnumcl=1:k
    tindxs = find(labels==tnumcl);
    sumLCPOp = sumLCPOp + sum(lcp(tindxs, tnumcl));
end

% Compute the IC scores
valIC.BIC = -2*llh + (C_phi_beta_k_bic * P_k);
valIC.AIC = -2*llh + (C_phi_beta_k_aic * P_k);
valIC.ICL = -2*llh + (log(N) * P_k) - (2*sumLCPOp);
valIC.EME = sumLCPOp;
valIC.LLHT = -2*llh;
valIC.COMT = (log(N) * P_k);
valIC.beta_min = -2*llh + (C_phi_beta_k_beta_min * P_k);
valIC.beta_max = -2*llh + (C_phi_beta_k_beta_max * P_k);
valIC.beta_range = -2*llh + (C_phi_beta_k_beta_range .* P_k);
valIC.beta_range_actVal = beta_range;
valIC.llh = llh;
end