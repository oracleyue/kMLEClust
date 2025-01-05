clear all; close all

addpath('../measures/')

% dataset name
name = 'Wafer'; K = 2;

% ground truth
fname = ['data/', name, '_labels.txt'];
labels_gt = readmatrix(fname, 'FileType', 'text', 'Delimiter', ' ')';

% use k-GAK (kernel kmeans)
methods = {'kMeansDTW', 'kShape', 'kGAK', 'kSC'};
for k = 1:length(methods)
    method = methods{k};

    fname = ['results/', name, '_', method, '_labels.txt'];
    labels = readmatrix(fname, 'FileType', 'text', 'Delimiter', ',');
    labels(~labels) = K;

    valRI = perfRI(labels, labels_gt, K);
    valNMI = perfNMI(labels, labels_gt, K);
    valARI = perfARI(labels, labels_gt, K);
    valNID = perfNID(labels, labels_gt, K);

    fprintf('%s:\n', method);
    fprintf('RI = %.4f\tNMI = %.4f\tARI = %.4f\t1-NID = %.4f\n\n',...
            valRI, valNMI, valARI, 1-valNID);
end
