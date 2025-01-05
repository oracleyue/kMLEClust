clear all; close all

% toolkit paths
addpath('../../toolkits/')
addpath('../../toolkits/k-SC/')
addpath('../../toolkits/u-shapelet/')
addpath('../../measures/')

% paths
datapath = '../data_expr1/';
respath  = '../result_expr1/';

% dataset spec.
mList = [2, 4, 8];
nPara = length(mList);
nExpr = 40;
K = 8;

% load dataset
load([datapath 'datasets_expr1.mat']);

% variables for results
labelKSC = cell(1, nPara);
precKSC = zeros(nExpr, nPara);  % Rand Index

%% benchmark k-SC
fprintf('Applying k-SC: \n');
for iPara = 1:nPara
    m = mList(iPara);
    fprintf('benchmark for m=%d: start\n', m)

    labelM = []; % each row is clustering result for a dataset
    prec = [];
    for iExpr = 1:nExpr
        data = dataExpr{iPara, iExpr};

        % init
        labels = labelExpr{iPara, iExpr};
        idxInit = zeros(1, K);
        for k = 1:K
            clust = find(labels == k);
            idxInit(k) = clust(1);
        end

        % clustering
        stime = tic;
        try
            groups = ksc_wrapper(data, K);
        catch ME
            switch ME.identifier
              case 'vecKARs:InvalidResponsibilities'
                groups = nan(N,1);
              otherwise
                rethrow(ME)
            end
        end

        % results
        accuracy = perfRI(groups, labels, K);
        etime = toc(stime);
        fprintf('... dataset #%d: cpu time %.4f, precision %.4f \n', ...
                iExpr, etime, accuracy);

        % saving result
        labelM = [labelM; groups'];
        prec = [prec; accuracy];
    end
    labelKSC{iPara} = labelM;
    precKSC(:, iPara) = prec;
end
% saving into .mat
save([respath, 'kSC_result.mat'], 'mList', 'labelKSC', 'precKSC');
