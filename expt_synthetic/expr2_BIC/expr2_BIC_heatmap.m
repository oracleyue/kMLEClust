clear all; close all

addpath('../../functions/')
addpath('../../measures/')

%% Model Simulation

% model specification
m = 4;            % m-variate VAR
p = 5;            % VAR(p)
T = 200;          % time length
K = 10;           % number of clusters
NperC = 20;       % number of time series per cluster
N = K * NperC;    % total number of time series

KList = 2:2:20;   % range of K for BIC
pList = 2:1:7;    % range of p for BIC

% data generation
rng(1)
[data, mdls, labels_gt] = simVARs(m, p, T, N, K);

% init
idxInit = zeros(1, K);
for k = 1:K
    clust = find(labels_gt == k);
    idxInit(k) = clust(1);
end
maxK = max(KList);
idxExtra = randi(N, 1, maxK-K);
idxInit = [idxInit, idxExtra];


%% Choose K by BIC in Clustering

fid_log = fopen('../result_expr2/BIC_heatmap.log', 'w');

% Vector k-ARs
scoreMatrix = nan(length(pList), length(KList));
precMatrix = nan(length(pList), length(KList));
for iK = 1:length(KList)
    Kval = KList(iK);

    for ip = 1:length(pList)
        pval = pList(ip);

        % clustering
        stimer = tic;
        [labels, mdl, loss] = kVARs(data, Kval, pval, 'init', idxInit(1:Kval));
        etime = toc(stimer);

        % precision
        valRI = perfRI(labels, labels_gt, K);
        precMatrix(ip, iK) = valRI;

        % BIC scoring
        score = calcBIC(data, mdl);
        scoreMatrix(ip, iK) = score;

        fprintf('(K=%d, p=%d) BIC score = %.2f, Prec = %.2f; %.4f sec. \n', ...
                Kval, pval, score, valRI, etime);
        fprintf(fid_log, '(K=%d, p=%d) BIC score = %.2f, Prec = %.2f; %.4f sec. \n', ...
                Kval, pval, score, valRI, etime);
    end
end
save('../result_expr2/BIC_heatmap.mat', 'KList', 'pList', 'scoreMatrix', 'precMatrix');


%% Visualization (heatmap of K & p)
fig_hl = figure(1);
h = heatmap(KList, pList, log(scoreMatrix));
h.XLabel = 'number of clusters (K)';
h.YLabel = 'model order (p)';
h.Title = 'BIC score in log';

% export in pdf
pos = [8.5417 10.5694 6.7778 3.5556];
set(fig_hl, 'Units','Inches', 'Position', pos);
set(fig_hl,'PaperPositionMode','Auto','PaperUnits','Inches',...
           'PaperSize',[pos(3), pos(4)]);
set(fig_hl, 'Renderer', 'Painters');  % enforce vector figure
print(fig_hl, 'BIC_heatmap.pdf', '-dpdf', '-r0')

%% end logging
fclose(fid_log)
