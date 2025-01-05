clear all; close all
rng(2)

% model specification
m = 3;            % m-variate VAR
p = 5;            % VAR(p)
T = 100;          % time length
K = 2;            % number of clusters
Nc = 20;          % number of ts per cluster
N = K * Nc;       % total number of time series

% data generation
addpath('../functions/')
[data, mdls, labels_gt] = simVARs(m, p, T, Nc, K, 'per_cluster');

% clustering
tic
[labels, centers] = kmeansVARs(data, K, p);
toc

% results
addpath('../measures/')
valRI  = perfRI(labels, labels_gt, K);
valARI = perfARI(labels, labels_gt, K);
valNMI = perfNMI(labels, labels_gt, K);
valNID = 1 - perfNID(labels, labels_gt, K);
fprintf('RI = %.2f\nARI = %.2f\nNMI = %.2f\n1-NID = %.2f\n',...
        valRI, valARI, valNMI, valNID)
