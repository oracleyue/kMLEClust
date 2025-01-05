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
tslist = [1 64 100 420];

fig_hl = figure;
y1 = data{tslist(1)};
y2 = data{tslist(2)};
y3 = data{tslist(3)};
y4 = data{tslist(4)};
cmap = colormap('lines');
for i = 1:m
    subplot(m/2, 2, i)
    plot(y1(i,:), '--^',  'color', cmap(2,:), 'linewidth', 1, 'markersize', 3.5)
    hold on
    plot(y2(i,:), '--v',  'color', cmap(2,:), 'linewidth', 1, 'markersize', 3.5)
    plot(y3(i,:), ':o',   'color', cmap(1,:), 'linewidth', 1, 'markersize', 3.5)
    plot(y3(i,:), ':*',   'color', cmap(1,:), 'linewidth', 1, 'markersize', 3.5)
    ylabel(['#' num2str(i) ' component series'])
    xlabel('Time')
    legend('#1, abnormal', '#64, abnormal', ...
           '#100, normal', '#420, normal')
    grid on
    hold off
end

% export in pdf
pos = [9.2222 9.2083 12.0278 6.5973];
set(fig_hl, 'Units','Inches', 'Position', pos);
set(fig_hl,'PaperPositionMode','Auto','PaperUnits','Inches',...
           'PaperSize',[pos(3), pos(4)]);
set(fig_hl, 'Renderer', 'Painters');  % enforce vector figure
print(fig_hl, ['ts_plots.pdf'], '-dpdf', '-r0')
