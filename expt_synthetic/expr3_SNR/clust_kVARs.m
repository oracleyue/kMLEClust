clear all; close all

addpath('../../functions/')
addpath('../../measures/')

% paths
datapath = '../data_expr3/';
respath  = '../result_expr3/';

% load dataset
load([datapath 'datasets_expr3.mat']);

% variables for results
labelKVARs = cell(1, nPara);
labelMVAR = cell(1, nPara);
precKVARs = zeros(nExpr, nPara);
precMVAR = zeros(nExpr, nPara);

%% benchmark k-VARs
fprintf('Applying k-VARs: \n');
for iPara = 1:nPara
    snr = SNRVec(iPara);
    fprintf('benchmark for SNR=%d dB: start\n', snr)

    labelM = []; % each row is clustering result for a dataset
    prec = [];
    for iExpr = 1:nExpr
        data = dataExpr{iPara, iExpr};

        % oracle init
        labels = labelExpr{iPara, iExpr};
        [~, idxInit] = unique(labels, 'first');

        % clustering
        simTimer = tic;
        [groups, mdl, loss] = kVARs(data, K, p, 'init', idxInit);
        etime = toc(simTimer);

        % results
        valPerf = perfNMI(groups, labels, K);
        fprintf('... dataset #%d: cpu time %.4f, precision %.4f%% \n', ...
                iExpr, etime, valPerf*100);

        % saving result
        labelM = [labelM; groups'];
        prec = [prec; valPerf];
    end
    labelKVARs{iPara} = labelM;
    precKVARs(:, iPara) = prec;
end
% saving into .mat
save([respath, 'kVARs_result.mat'], 'SNRVec', 'labelKVARs', 'precKVARs');
fprintf('\n')
