clear all; close all

% toolkit paths
addpath('../toolkits/')
addpath('../toolkits/k-SC/')

% paths
respath = 'results/';

% load dataset
load('./data/Wafer.mat');

rawdata = mts.test;
labels = mts.testlabels;
K = 2;

% data reformating
T = 104;
data = {};
for n = 1:length(rawdata)
    item = rawdata{n};
    item = item(:, 1:T);
    data = [data, item];
end

% clustering
groups = ksc_wrapper(data, K);

%% saving results
fid = fopen([respath 'Wafer_kSC_labels.txt'], 'w');
for i = 1:length(groups)
    fprintf(fid, '%d\n', groups(i));
end
fclose(fid);
