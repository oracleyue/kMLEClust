clear all; close all

% toolkit paths
addpath('../../toolkits/')
addpath('../../toolkits/k-SC/')
addpath('../../toolkits/u-shapelet/')
addpath('../../measures/')

% paths
datapath = '../data_expr3/';
respath  = '../result_expr3/';

% load dataset
load([datapath 'datasets_expr3.mat']);

% variables for results
labelKSC = cell(1, nPara);
precKSC = zeros(nExpr, nPara);

%% benchmark k-SC
fprintf('Applying k-SC: \n');
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
        groups = ksc_wrapper(data, K);
        etime = toc(simTimer);

        % results
        valPerf = perfNMI(groups, labels, K);
        fprintf('... dataset #%d: cpu time %.4f, precision %.4f%% \n', ...
                iExpr, etime, valPerf*100);

        % saving result
        labelM = [labelM; groups'];
        prec = [prec; valPerf];
    end
    labelKSC{iPara} = labelM;
    precKSC(:, iPara) = prec;
end
%% saving into .mat
save([respath, 'kSC_result.mat'], 'SNRVec', 'labelKSC', 'precKSC');
fprintf('\n')
