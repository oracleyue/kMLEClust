clear all; close all

addpath('../measures/')

% load results
respath = './results/';
load([respath 'wafer_results_kVARs.mat'])

% compute perfmance indexes
N = length(cGroups);
valRI = zeros(1, N);
valARI = zeros(1, N);
valNMI = zeros(1, N);
valNID = zeros(1, N);

for n = 1:N
    groups = cGroups{n};
    valRI(n)  = perfRI (groups, labels, K);
    valARI(n) = perfARI(groups, labels, K);
    valNMI(n) = perfNMI(groups, labels, K);
    valNID(n) = 1 - perfNID(groups, labels, K);
end

% quartile summary
stats = quantile(valRI,[0.25 0.50 0.75]);
fprintf('RI:\t %f, %f, %f\n', stats(1), stats(2), stats(3));
stats = quantile(valNMI,[0.25 0.50 0.75]);
fprintf('NMI:\t %f, %f, %f\n', stats(1), stats(2), stats(3));
stats = quantile(valARI,[0.25 0.50 0.75]);
fprintf('ARI:\t %f, %f, %f\n', stats(1), stats(2), stats(3));
stats = quantile(valNID,[0.25 0.50 0.75]);
fprintf('1-NID:\t %f, %f, %f\n', stats(1), stats(2), stats(3));

% boxplots
fhl = figure(1);
boxplot([valRI', valNMI', valARI', valNID'], 'Notch','off',...
        'Labels',{'RI', 'NMI', 'ARI', '1-NID'});
grid on

pos = [8.5417 12.4861 4.8889 1.6389];
set(fhl, 'units', 'inches', 'position', pos);
set(fhl,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fhl, [respath 'Wafer_kVARs_boxplot.pdf'], '-dpdf', '-r0')
