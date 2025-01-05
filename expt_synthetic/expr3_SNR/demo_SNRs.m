clear all; close all
rng(2)

% paths
addpath('../../functions/')
addpath('../../measures/')

% model specification
p = 2;            % VAR(p)
m = 3;            % variate dim
T = 200;          % time length
K = 4;            % number of clusters
Nc = 40;          % number of ts per cluster
N = K * Nc;       % total number of time series

% data generation
SNRList = [-2, 0, 2, 4, 8];
[dataList, mdlList_gt, labelList_gt] = simVARs_SNRs(m, p, T, SNRList, N, K);

% choose data with specific SNR
snr = -2;
idx = find(SNRList == snr);
data = dataList{idx};
labels_gt = labelList_gt{idx};

% random init
% [labels, mdl, loss] = kVARs(data, K, p);
% oracle init
[~, idxInit] = unique(labels_gt, 'first');

% clustering
tic
[labels, mdl, loss] = kVARs(data, K, p, 'init', idxInit);
toc

% results
valRI  = perfRI( labels, labels_gt, K);
valARI = perfARI(labels, labels_gt, K);
valNMI = perfNMI(labels, labels_gt, K);
valNID = perfNID(labels, labels_gt, K);
disp('k-VARs:')
fprintf('RI = %.2f\nARI = %.2f\nNMI_root = %.2f\nNMI_max = %.2f\n',...
        valRI, valARI, valNMI, 1-valNID)
