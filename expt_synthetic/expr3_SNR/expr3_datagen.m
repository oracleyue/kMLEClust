clear all; close all

% search path
addpath('../../functions/')
% data path
datapath = '../data_expr3/';

% seeding
rng(2)

% model specification
m = 3;                      % m-variate VAR(p)
p = 5;                      % model order or lags
SNRVec = [0, 2, 4, 8, 10];  % signal-noise-ratio (dB)
T = 80;                     % time length
K = 8;                      % number of clusters
Nc = 30;                    % number of ts per cluster
N = K * Nc;                 % total number of ts
nExpr = 40;                 % number of experiments

% data generation
nPara = length(SNRVec);
dataExpr  = cell(nPara, nExpr);
mdlExpr   = cell(nPara, nExpr);
labelExpr = cell(nPara, nExpr);
SNRExpr = cell(1, nExpr);

% logging
diary modelGeneration.log

%% generating data
fprintf('Data generating:\n\n')
simTimer = tic;
for j = 1:nExpr
    [dataList, mdlList, labelList, snrList] = ...
        simVARs_SNRs(m, p, T, SNRVec, Nc, K, 'per_cluster');

    for k = 1:nPara
        dataExpr{k,j}  = dataList{k};
        mdlExpr{k,j}   = mdlList{k};
        labelExpr{k,j} = labelList{k};
    end
    SNRExpr{j} = snrList;
    fprintf('... #%d dataset (with %d SNRs): done\n\n\n', j, nPara);
end
etime = toc(simTimer);
fprintf('model simulation done: %d models, %.2f sec.\n', nPara*nExpr, etime);
save([datapath 'datasets_expr5.mat'])

% logging off
diary off

%% export data in txt for tslearn Python package
for iPara = 1:nPara
    snr = SNRVec(iPara);

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
    for iExpr = 1:nExpr
        data = dataExpr{iPara, iExpr};

        fname = [datapath sprintf('SNR%d_dataset%d.txt', snr, iExpr)];
        fid = fopen(fname, 'w');
        for n = 1:N
            ts = data{n}';
            fprintf(fid, formatSpec, ts);
            if n < N
                fprintf(fid, '\n');
            end
        end
        fclose(fid);
    end
end

%% export ground truth labels in txt
for iPara = 1:nPara
    snr = SNRVec(iPara);

    % formating spec.
    formatSpec = '';
    for j = 1:N
        if j == 1
            formatSpec = [formatSpec '%d'];
        else
            formatSpec = [formatSpec ' %d'];
        end
    end

    % writing to file
    for iExpr = 1:nExpr
        data = labelExpr{iPara, iExpr};

        fname = [datapath sprintf('SNR%d_dataset%d_labels.txt', snr, iExpr)];
        fid = fopen(fname, 'w');
        fprintf(fid, formatSpec, data);
        fclose(fid);
    end
end
