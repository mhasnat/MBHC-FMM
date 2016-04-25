% This file performs the following tasks:
% 1. Load image normals and
% 2. Cluster normals using MBHC-FMM method (for depth image analysis)

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

% Load depth image
load('depImg.mat');
imgNormals = depImg.imgNormals;
[r, c, d] = size(imgNormals);

% Display image normals
subplot(121);
imshow(imgNormals, []);
title('Image normals');

% reshape as (N x 3)
imgNormals = reshape(imgNormals, r*c, d);

%% Part 2 (data clustering - using MBHC-FMM method)
% Perform clustering data
kMax = 15; % desired maximum possible number of clusters

% get model parameters with k_max
params = bd_vmfmm(imgNormals, kMax);
params = annihilateComp(params); % eliminate empty clusters

% get model parameters with k_max
[~, allIC, allClust] = get_HC_IC_BD(imgNormals, params);

% Model selection
y = allIC.BIC;
x = 1: length(y);

% Select the number of clusters automatically (s.t. k_max number of
% clusters)
[numComp,~, ~] = wplr(x,y,[1 30]); % WLR-30

% Get the final clustering results
finalClust = allClust(:, numComp);

% Display data in sphere with label
labelImg = reshape(finalClust, r, c);
subplot(122);
imshow(label2rgb(labelImg), []);
title('Depth Image Regions');