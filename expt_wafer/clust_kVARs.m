clear all; close all
rng(4)

% addpath
addpath('../functions/')
addpath('../measures/')

% load data
load('./data/Wafer.mat')
data = mts.test;
labels = mts.testlabels;

% model specification
m = 6;            % m-variate VAR
p = 10;           % VAR(p); 10
T = 104;          % time length
K = 2;            % number of clusters
N = 896;          % total number of time series

nRuns = 20;     % the number of results for stat. plot

% saving results
cGroups = {};
cLld = {};
initPairs = zeros(2, nRuns);

% multi-running
ns = 0; nr = 0;
while nr < nRuns
    ns = ns + 1;
    fprintf('\nRunning with random init. No.%d:\n', ns);

    % clustering
    idxInit = [randi(95), randi([96,N])];
    % idxInit = randi(N, [1 2]);
    [groups, mdl, loss] = kVARs(data, K, p, 'init', idxInit,...
                                'maxIter', 200, 'tol', 1e-10);

    % display results
    valRI  = perfRI (groups, labels, K);
    valARI = perfARI(groups, labels, K);
    valNMI = perfNMI(groups, labels, K);
    valNID = perfNID(groups, labels, K);
    fprintf('RI = %.4f\tARI = %.4f\tNMI = %.4f\t1-NID = %.4f\n',...
            valRI, valARI, valNMI, 1-valNID)
    fprintf('loss (log) = %.6f\n', log(loss))

    % thresholding likelihood
    % Note: the threshold value can be chosen by firstly dozens of random runs.
    lldThreshold = 15.57;
    if log(loss) > lldThreshold
        display('hit!')
        nr = nr + 1;
        cGroups = [cGroups, groups];
        cLld = [cLld, log(loss)];
        initPairs(:, nr) = idxInit';
    end
end

%% saving results
save('./results/wafer_results_kVARs.mat', 'labels', 'K', 'cGroups', 'cLld', 'initPairs');
