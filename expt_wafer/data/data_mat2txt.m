clear all; close all

% load data
load('Wafer.mat')
data_summary  % print info of dataset

% data = mts.train;
% labels = mts.trainlabels;
data = mts.test;
labels = mts.testlabels;
% mtsC = combineData(mts);
% data = mtsC.data;
% labels = mtsC.labels;

% data specification
m = 6;            % m-variate VAR
p = 8;            % VAR(p)
T = 104;          % time length
K = 2;            % number of clusters
N = 896;          % total number of time series

% reformat and export data in txt for tslearn Python package
% formating spec.
formatSpec = '';
for i = 1:m
    for j = 1:T
        if j == 1
            formatSpec = [formatSpec '%f'];
        else
            formatSpec = [formatSpec ' %f'];
        end
    end
    if i < m
        formatSpec = [formatSpec '|'];
    end
end

% writing to file
fname = 'Wafer_data.txt';
fid = fopen(fname, 'w');
for n = 1:N
    ts = data{n}';
    ts = ts(1:T, 1:m);
    fprintf(fid, formatSpec, ts);
    if n < N
        fprintf(fid, '\n');
    end
end
fclose(fid);

fname = 'Wafer_labels.txt';
writematrix(labels', fname, 'FileType', 'text', 'Delimiter', ' ');
