function params = bd_vmfmm(vectors, k, clust)
% Performs clustering using vMF mixture model. 
% See Sect. 3 and 4 of ref [1] or Sect. 2 of ref [2]

% INPUT:
% vectors: feature vectors (N x 3)
% k      : Number of clusters
% clust (optional): Initial labels (cluster assignments) for the data

% OUTPUT
% params: different types of information related to clustering.
% *note* If the the cluster labels are desired, then one can directly get them
% from 'params.label' and ignore the rest

% Reference:
% [1] Hasnat et al., Model-based hierarchical clustering with Bregman 
% divergences and Fishers mixture model: application to depth image analysis. 
% Statistics and Computing, 1-20, 2015.
%
% [2] Hasnat et al., Hierarchical 3-D von Mises-Fisher Mixture Model, ICML-WDDL, 2013.
% 

% Author: Md Abul HASNAT

desNumComp = k;
[D,V] = size(vectors);
numOfDataSample = D;
dim   = V;

if(nargin<3)
    clust = kmeanspp(vectors',k);
end

%% Bragman soft clustering
MAX_ITERATIONS = 200;

% Initialization
% Compute the weight and expectation parameter from the available
% information for each cluster
i = 1;
for j=1:desNumComp
    numOfDataPointToThisCluster = length(find(clust==j));
    
    alpha(i) = numOfDataPointToThisCluster / numOfDataSample; % weight
    
    % Expectation parameter computation
    indx = find(clust==j);
    dataClust = vectors(indx, :);
    
    sufStat(i, :) = sum(dataClust);
    
    eta(i, :) = sufStat(i, :) / length(dataClust);
    normEta(i) = sqrt(eta(i, :) * eta(i, :)');
    normTheta(i) = getThetaFromEta(normEta(i));
    
    % Compute R(normTheta)
    R_norm_theta(i) = ((1/tanh(normTheta(i))) - (1/normTheta(i))) / normTheta(i);
    theta_cl(i, :) = eta(i, :) ./ R_norm_theta(i); % natural parameter
    i = i+1;
end
clear i;

% Start the EM algorithm
logLikelihoodNew = 100;
logLikelihoodOld = 0;
logLikelihoodThreshold = 0.001;

iterations = 1;
diffLLH = abs(logLikelihoodOld-logLikelihoodNew);
while(iterations<MAX_ITERATIONS && diffLLH>logLikelihoodThreshold)
    logLikelihoodOld = logLikelihoodNew;
    
    % Expectation Step
    for j=1:k
        % Compute Bragman divergence
        Log_Normalizing_Function(j) = log((4*pi*sinh(normTheta(j))) / normTheta(j));
        Dual_Log_Normalizing_Function(j) = (eta(j, :) * theta_cl(j, :)') - Log_Normalizing_Function(j);
        Grad_Dual_Log_Normalizing_Function(j,:) = theta_cl(j, :);
        
        sufStat_minus_expectParam(:, :, j) = bsxfun(@minus, vectors , eta(j, :));
        innerProdTerm(:,j) = sufStat_minus_expectParam(:, :, j) * Grad_Dual_Log_Normalizing_Function(j,:)';
        
        % Compute divergence
        divergence(:,j) = -(Dual_Log_Normalizing_Function(j) + innerProdTerm(:,j));
        expTerm(:,j) = exp(-divergence(:,j));
        probTerm(:, j) = alpha(j) * expTerm(:,j);
    end
    
    probTerm = bsxfun(@rdivide, probTerm, sum(probTerm, 2));
    [~, clust] = max(probTerm,[], 2);
    
    % Maximization Step
    sumProbTerm = sum(probTerm);
    
    % Update weight
    alpha = sumProbTerm ./ size(probTerm,1);
    
    % Update parameter
    for j=1:k
        eta(j, :) = sum(bsxfun(@times, probTerm(:,j) , vectors)) ./ sumProbTerm(j);
        
        % Convert to source parameters (mu, kappa) from expectation parameter
        % (eta)
        normEta(j) = sqrt(eta(j, :) * eta(j, :)');
        normTheta(j) = getThetaFromEta(normEta(j));
        
        % Compute R(normTheta)
        R_norm_theta(j) = ((1/tanh(normTheta(j))) - (1/normTheta(j))) / normTheta(j);
        theta_cl(j, :) = eta(j, :) ./ R_norm_theta(j);
        
        update_kappa(j) = normTheta(j);
        update_mu(j, :) = theta_cl(j, :) ./ normTheta(j);
    end
    
    % Compute new log likelihood
    logWeight   = log(alpha);
    logNormTerm = log(update_kappa) - log(4*pi*sinh(update_kappa));
    logExpTerm  = bsxfun(@times, update_kappa,  (update_mu * vectors')');
    
    logClassCondLiklihood = bsxfun(@plus, logWeight + logNormTerm , logExpTerm);
    ClassCondLiklihood    = exp(logClassCondLiklihood);
    logMixtureLikelihoodPerComponent = log(sum(ClassCondLiklihood,2));
    logLikelihoodNew = sum(logMixtureLikelihoodPerComponent);
    diffLLH = abs(logLikelihoodOld-logLikelihoodNew);
    
    LLH(iterations) = -logLikelihoodNew;
    
    iterations = iterations+1;
end

%% Save Parameters
params.expectation = eta;
params.natural = theta_cl;
params.source.kappa = update_kappa;
params.source.mu = update_mu;
params.weight = alpha;
params.label = clust;
params.cp = probTerm;
