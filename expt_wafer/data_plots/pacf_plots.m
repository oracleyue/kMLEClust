clear all; close all
rng(2)

% load data
load('../data/Wafer.mat')

addpath('../data/')
data_summary  % print info of dataset

% data = mts.train;
% labels = mts.trainlabels;
data = mts.test;
labels = mts.testlabels;
% data = [mts.train mts.test];
% labels = [mts.trainlabels; mts.testlabels];

% model specification
m = 6;            % m-variate VAR
p = 8;            % VAR(p)
T = 104;          % time length
K = 2;            % number of clusters
N = 896;          % total number of time series

%% plot ACF
tslist = [1 64 100 420]; % 1 64 (abnormal); 100 420 (normal)

fig_hl = figure;
for idx = 1:length(tslist)
    y = data{tslist(idx)};
    for j = 1:m
        subplot(4, m, (idx-1)*m+j)
        parcorr(y(j,:), 'NumLags',20, 'NumSTD',2)
        ylabel(['#' num2str(j) ' PACF'])
        xlabel('Lag')
        title('')
    end
end

% export in pdf
% set(fig_hl, 'Units','Inches');
pos = [7.6250 3.9028 19.0694 10];
set(fig_hl, 'Units','Inches', 'Position', pos);
set(fig_hl,'PaperPositionMode','Auto','PaperUnits','Inches',...
           'PaperSize',[pos(3), pos(4)]);
set(fig_hl, 'Renderer', 'Painters');  % enforce vector figure
print(fig_hl, ['pacf_plots.pdf'], '-dpdf', '-r0')
