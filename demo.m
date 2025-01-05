clear all; close all
rng(2)

addpath('./functions/')
addpath('./measures/')

% model specification
m = 1;            % m-variate VAR
p = 5;            % VAR(p)
T = 100;          % time length
K = 2;            % number of clusters
Nc = 20;          % number of ts per cluster
N = K * Nc;       % total number of time series

% data generation
[data, mdls, labels_gt] = simVARs(m, p, T, Nc, K, 'per_cluster');

% clustering
tic
index = 2;
switch index
    case 1
        % random initialization
        [labels, mdl, loss] = kVARs(data, K, p, 'initType', 'random');
    case 2
        % multiple random initialization (automatically set: 'initType', 'random')
        [labels, mdl, loss] = kVARs(data, K, p, 'init', 5);
    case 3
        % initialization via kmeans
        [labels, mdl, loss] = kVARs(data, K, p, 'initType', 'kmeans');
    case 4
        % user-specified initialization (automatically set: 'initType', 'user')
        idxInit = [1, 30];
        [labels, mdl, loss] = kVARs(data, K, p, 'init', idxInit);
end
toc

% results
valRI  = perfRI(labels, labels_gt, K);
valARI = perfARI(labels, labels_gt, K);
valNMI = perfNMI(labels, labels_gt, K);
valNID = 1 - perfNID(labels, labels_gt, K);
fprintf('RI = %.2f\nARI = %.2f\nNMI = %.2f\n1-NID = %.2f\n',...
        valRI, valARI, valNMI, valNID)
