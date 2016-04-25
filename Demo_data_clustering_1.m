% This file performs the following tasks:
% 1. Load pre-generated samples/data and
% 2. Cluster the data using MBHC-FMM method

% For the GUI version, see mbc_vmfmm.m

% Reference:
% [1] Hasnat et al., Model-based hierarchical clustering with Bregman 
% divergences and Fishers mixture model: application to depth image analysis. 
% Statistics and Computing, 1-20, 2015.
%
% [2] Hasnat et al., Hierarchical 3-D von Mises-Fisher Mixture Model, ICML-WDDL, 2013.
% 

% Author: Md Abul HASNAT

clc; clear all; close all;

%% Part 1 (Load data and ground truth labels)

% Load samples and labels
FileName = 'vmfmm_samp.mat';
load(FileName);

% Display data in sphere with label
spread(vmfSample', labels);
title('Samples distributed w.r.t. vMF distribution');

%% Part 2 (data clustering - using MBHC-FMM method)
% Perform clustering data
kMax = 15; % desired maximum possible number of clusters

% get model parameters with k_max
params = bd_vmfmm(vmfSample, kMax);
params = annihilateComp(params); % eliminate empty clusters

% get model parameters with k_max
[~, allIC, allClust] = get_HC_IC_BD(vmfSample, params);

% Model selection
y = allIC.BIC;
x = 1: length(y);

% Select the number of clusters automatically (s.t. k_max number of
% clusters)
[numComp,~, ~] = wplr(x,y,[1 30]); % WLR-30

% Get the final clustering results
finalClust = allClust(:, numComp);

% Display data in sphere with label
figure;
spread(vmfSample', finalClust);
title('Samples clustered with MBHC-FMM method');