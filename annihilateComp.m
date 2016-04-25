function params = annihilateComp(params)
% Performs component annihilation/deletion. That means, removes invalid
% clusters/components of a vMF mixture model

% INPUT:
% params: parameters of a vMF mixture model

% OUTPUT
% params: parameters of a vMF mixture model after removal of invalid clusters

% Reference:
% [1] Hasnat et al., Model-based hierarchical clustering with Bregman 
% divergences and Fishers mixture model: application to depth image analysis. 
% Statistics and Computing, 1-20, 2015.
%
% [2] Hasnat et al., Hierarchical 3-D von Mises-Fisher Mixture Model, ICML-WDDL, 2013.
% 

% Author: Md Abul HASNAT

wt = params.weight;

indx = [];
for i=1:length(wt)
    if(wt(i)<0.0005)
        indx = [indx; i];
    end
end

if(~isempty(indx))
    params.cp(:,indx) = [];
    params.expectation(indx,:) = [];
    params.natural(indx,:) = [];
    params.weight(indx) = [];
    
    params.source.kappa(indx) = [];
    params.source.mu(indx,:) = [];
end