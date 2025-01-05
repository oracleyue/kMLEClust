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
tsIdx = 64;    % 1 64; 100 420 (normal)

fig_hl = figure;
y = data{tsIdx};
for i = 1:m
    for j = 1:m
        subplot(m,m, (i-1)*m+j)
        if i == j
            autocorr(y(i,:), 'NumLags',20, 'NumSTD',2)
            ylabel('ACF')
        else
            crosscorr(y(i,:), y(j,:), 'NumLags',20, 'NumSTD',2)
            ylabel('XCF')
        end

        xlabel('Lag')
        title('')
    end
end

% export in pdf
pos = [5.2361 3.8194 21.4167 11.3889];
set(fig_hl, 'Units','Inches', 'Position', pos);
set(fig_hl,'PaperPositionMode','Auto','PaperUnits','Inches',...
           'PaperSize',[pos(3), pos(4)]);
set(fig_hl, 'Renderer', 'Painters');  % enforce vector figure
print(fig_hl, ['acf_plot_' num2str(tsIdx) '.pdf'], '-dpdf', '-r0')
