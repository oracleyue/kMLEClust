clear all; close all
rng(2)

% toolkit paths
addpath('../functions/')
addpath('./k-SC/')
addpath('./u-shapelet/')

% model specification
m = 9;            % m-variate VAR
p = 5;            % VAR(p)
T = 100;          % time length
K = 3;            % number of clusters
Nc = 10;          % number of ts per cluster
N = K * Nc;       % total number of time series

% data generation
[data, ~, labels] = simVARs(m, p, T, Nc, K, 'per_cluster');

% clustering
tic
groups = kSC_wrapper(data, K);
toc

% results
valRI = perfRI(groups, labels, K);
valNMI = perfNMI(groups, labels, K);
fprintf('RI: %.4f    NMI: %.4f\n', valRI, valNMI);
