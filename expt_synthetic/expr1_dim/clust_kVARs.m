clear all; close all

addpath('../../functions/')
addpath('../../measures/')

% paths
datapath = '../data_expr1/';
respath  = '../result_expr1/';

% dataset spec.
mList = [2, 4, 8];
nPara = length(mList);
nExpr = 40;
K = 8;

% method spec.
type = 'rand';
% type = 'oracle';

% load dataset
load([datapath 'datasets_expr1.mat']);

% variables for results
labelKVARs = cell(1, nPara);
labelMVAR = cell(1, nPara);
precKVARs = zeros(nExpr, nPara);
precMVAR = zeros(nExpr, nPara);

%% benchmark k-VARs
fprintf('Applying k-VARs: \n');
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
        switch type
            case 'rand'
                % use random initialization
                try
                    [groups, mdl, loss] = kVARs(data, K, p, 'initType', 'random', 'init', 10);
                catch ME
                    switch ME.identifier
                        case 'kVARs:BadInitialization'
                            [groups, mdl, loss] = kVARs(data, K, p, 'init', idxInit);
                        otherwise
                            rethrow(ME)
                    end
                end
            case 'oracle'
                [groups, mdl, loss] = kVARs(data, K, p, 'init', idxInit);
            otherwise
                error('Not supported!')
        end

        % results
        accuracy = perfNMI(groups, labels, K);
        etime = toc(stime);
        fprintf('... dataset #%d: cpu time %.4f, accuracy %.4f%% \n', ...
                iExpr, etime, accuracy*100);

        % saving result
        labelM = [labelM; groups'];
        prec = [prec; accuracy];
    end
    labelKVARs{iPara} = labelM;
    precKVARs(:, iPara) = prec;
end

% saving into .mat
switch type
    case 'rand'
        save([respath, 'kVARs_result_rnd.mat'], 'mList', 'labelKVARs', 'precKVARs');
    case 'oracle'
        save([respath, 'kVARs_result.mat'], 'mList', 'labelKVARs', 'precKVARs');
end
