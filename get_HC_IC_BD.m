function [Clusters_Params, allIC_BD, allClust_BD] = get_HC_IC_BD(vectors, params)
% Construct a hierarchy of vMF mixture models by exploiting Bregman
% Divergence and compute the scores of different information criteria
% See Section 4 of [1] for details

% INPUT:
% vectors: feature vectors (N x 3)
% params: parameters of a vMF mixture model (with k_max
% components/clusters)

% OUTPUT
% Clusters_Params: Parameters of the vMF mixture models with different
% value of k (i.e., k = k_max-1, ..., 1)
% allIC_BD: scores of different information criteria
% allClust_BD: corresponding cluster labels for different
% value of k (i.e., k = k_max-1, ..., 1)

% Reference:
% [1] Hasnat et al., Model-based hierarchical clustering with Bregman 
% divergences and Fishers mixture model: application to depth image analysis. 
% Statistics and Computing, 1-20, 2015.

% Author: Md Abul HASNAT

% Get the parameters from the structure
eta = params.expectation;
theta_cl = params.natural;
weight = params.weight';
k = length(weight);
cp = params.cp;
label = params.label;

%% Perform hierarchical clustering
numberOfCluster = k;

LNF = zeros(numberOfCluster,1);
DLNF = zeros(numberOfCluster,1);

for j = 1:numberOfCluster
    normEta(j) = sqrt(eta(j, :) * eta(j, :)');
    normTheta(j) = getThetaFromEta(normEta(j));
    
    LNF(j) = log((4*pi*sinh(normTheta(j))) / normTheta(j));
    DLNF(j) = (eta(j, :) * theta_cl(j, :)') - LNF(j);
end
clear j;

% Compute distance among clusters
indx = 1;
Div_G = zeros(1, (numberOfCluster*(numberOfCluster-1))/2);
distMat = zeros(size(Div_G));

for i = 1 : numberOfCluster-1
    for j = i+1 : numberOfCluster
        % Left sided distance
        Div_G(indx) = DLNF(i) - DLNF(j) - ( (eta(i,:) - eta(j,:)) * theta_cl(j,:)');
        distMat(indx) = weight(i) * weight(j) * Div_G(indx);
        indx = indx+1;
    end
end
clear i j;

% Apply hierarchical clustering (use matlab function 'linkage')
Z = linkage(distMat,'average'); % 'average distance'
Clusters_Params = cell(1, numberOfCluster);

%% Perform clustering
for numClust = numberOfCluster : -1 : 1
    ucp = zeros(size(vectors,1), numClust);
    
    % Do initial clustering
    T = cluster(Z,'maxclust',numClust);
    isConsidered = zeros(size(T));
    
    % Allocate updated cluster parameter vectors
    Up_eta = zeros(numClust, size(eta,2));
    Up_theta_cl = zeros(numClust, size(theta_cl,2));
    Up_weight = zeros(numClust, 1);
    update_mu = zeros(numClust, size(eta,2));
    update_kappa = zeros(1,numClust);
    
    subset_cl = cell(1, numClust);
    indx=1;
   
    numClMerge = 0;
        
    % Find the merged clusters and Update cluster centroid (usign Bregman centroid)
    for i=1:numberOfCluster
        if(isConsidered(i)), continue; end

        cl_label = T(i);
        % get the associated subsets of original cluster
        subset_cl{indx} = find(T==cl_label);
        isConsidered(subset_cl{indx}) = 1;
        
        % Merge the subset and update centroid and weight
        if(length(subset_cl{indx}) > 1)
            % Left sided Bregman centroid computation with expectation
            % parameters
            Up_weight(indx) = sum(weight(subset_cl{indx}));
            Up_eta(indx,:) = sum(bsxfun(@times, weight(subset_cl{indx}), eta(subset_cl{indx},:))) ./ Up_weight(indx);
            
            % Convert to source parameters (mu, kappa) from expectation parameter
            % (eta)
            normUp_eta(indx) = sqrt(Up_eta(indx, :) * Up_eta(indx, :)');
            normUp_theta(indx) = getThetaFromEta(normUp_eta(indx));
            
            % Compute R(normTheta)
            R_norm_Up_theta(indx) = ((1/tanh(normUp_theta(indx))) - (1/normUp_theta(indx))) / normUp_theta(indx);
            Up_theta_cl(indx, :) = Up_eta(indx, :) ./ R_norm_Up_theta(indx);
            
            update_kappa(indx) = normUp_theta(indx);
            update_mu(indx, :) = Up_theta_cl(indx, :) ./ normUp_theta(indx);
            
            ucp(:, indx) = sum(cp(:, subset_cl{indx}),2);
            
            for ijk = 1:length(subset_cl{indx})
                numClMerge = numClMerge + length(find(label == subset_cl{indx}(ijk)));
            end
       else
            % no update
            Up_weight(indx) = weight(subset_cl{indx});
            Up_eta(indx,:) = eta(subset_cl{indx},:);
            
            % Convert to source parameters (mu, kappa) from expectation parameter
            % (eta)
            normUp_eta(indx) = sqrt(Up_eta(indx, :) * Up_eta(indx, :)');
            normUp_theta(indx) = getThetaFromEta(normUp_eta(indx));
            
            % Compute R(normTheta)
            R_norm_Up_theta(indx) = ((1/tanh(normUp_theta(indx))) - (1/normUp_theta(indx))) / normUp_theta(indx);
            Up_theta_cl(indx, :) = Up_eta(indx, :) ./ R_norm_Up_theta(indx);
            
            update_kappa(indx) = normUp_theta(indx);
            update_mu(indx, :) = Up_theta_cl(indx, :) ./ normUp_theta(indx);
            
            ucp(:, indx) = cp(:, subset_cl{indx});
        end
        
        indx = indx + 1;
    end
    
    Clusters_Params{numClust}.src_params.mu = update_mu;
    Clusters_Params{numClust}.src_params.kappa = update_kappa;
    Clusters_Params{numClust}.weight = Up_weight;
    Clusters_Params{numClust}.exp_params = Up_eta;
    Clusters_Params{numClust}.nat_params = Up_theta_cl;
    
    
    %% Compute Information Criterion values
    clear prms
    prms.mu = update_mu;
    prms.kappa = update_kappa;
    prms.weight = Up_weight';
    
    % --> Get the class number for each component
    clear sufStat_minus_expectParam Log_Normalizing_Function gradDLNF innerProdTerm divergence
    
    for j=1:numClust
        normEta2(j) = sqrt(Up_eta(j, :) * Up_eta(j, :)');
        normTheta2(j) = getThetaFromEta(normEta2(j));
        R_norm_theta2(j) = ((1/tanh(normTheta2(j))) - (1/normTheta2(j))) / normTheta2(j);
        theta_cl2(j, :) = Up_eta(j, :) ./ R_norm_theta2(j);
        
        Log_Normalizing_Function(j) = log((4*pi*sinh(normTheta2(j))) / normTheta2(j));
        Dual_Log_Normalizing_Function(j) = (Up_eta(j, :) * theta_cl2(j, :)') - Log_Normalizing_Function(j);
        gradDLNF(j,:) = theta_cl2(j, :);
        
        sufStat_minus_expectParam(:, :, j) = bsxfun(@minus, vectors , Up_eta(j, :));
        innerProdTerm(:,j) = sufStat_minus_expectParam(:, :, j) * gradDLNF(j,:)';
        divergence(:,j) = -(Dual_Log_Normalizing_Function(j) + innerProdTerm(:,j));
    end
       
    % Hard clustering
    [~, allClust_BD(:, numClust)] = min(divergence,[], 2);
    
    valIC2 = getICvalues_phi_beta_vmfmm(vectors, prms, allClust_BD(:, numClust), ucp);
    
    allIC_BD.BIC(numClust) = valIC2.BIC;
    allIC_BD.ICL(numClust) = valIC2.ICL;
    allIC_BD.AIC(numClust) = valIC2.AIC;
    allIC_BD.Beta_min(numClust) = valIC2.beta_min;
    allIC_BD.Beta_range(numClust,:) = valIC2.beta_range;
    allIC_BD.Beta_range_actVal(numClust,:) = valIC2.beta_range_actVal;
    
    allIC_BD.tE(numClust) = -sum(sum(ucp .* log2(ucp)));
    allIC_BD.numMerge(numClust) = numClMerge;    
end