clear all; close all

% search path
addpath('../../functions/')

% data path
datapath = '../data_expr1/';

% seeding
rng(4)

% model specification
mVec = [2 4 8];         % m-variate VAR(p)
p = 5;                  % model order or lags
T = 80;                 % time length
K = 8;                  % number of clusters
Nc = 30;                % number of ts per cluster
N = K * Nc;             % total number of ts
nExpr = 40;             % number of experiments

% data generation
nPara = length(mVec);
dataExpr  = cell(nPara, nExpr);
mdlExpr   = cell(nPara, nExpr);
labelExpr = cell(nPara, nExpr);

%% generating data
fprintf('Data generating:\n')
simTimer = tic;
for k = 1:nPara
    m = mVec(k);
    for j = 1:nExpr
        [data, mdl, labels] = simVARs(m, p, T, Nc, K, 'per_cluster');

        dataExpr{k,j}  = data;
        mdlExpr{k,j}   = mdl;
        labelExpr{k,j} = labels;

        fprintf('... #%d dataset (m=%d): done\n', j, m);
    end
end
etime = toc(simTimer);
fprintf('model simulation done: %d models, %.2f sec.\n', nPara*nExpr, etime);
save([datapath 'datasets_expr1.mat'])

%% export data in txt for tslearn Python package
for iPara = 1:nPara
    m = mVec(iPara);

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

        fname = [datapath sprintf('m%d_dataset%d.txt', m, iExpr)];
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
    m = mVec(iPara);

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

        fname = [datapath sprintf('m%d_dataset%d_labels.txt', m, iExpr)];
        fid = fopen(fname, 'w');
        fprintf(fid, formatSpec, data);
        fclose(fid);
    end
end
