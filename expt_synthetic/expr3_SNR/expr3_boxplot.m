clear all; close all

addpath('../../measures/')

% Retrieve dataset specifications
load('../data_expr3/datasets_expr3.mat')

% Choose performance index
% perfIndex = 'RI';
% perfIndex = 'NMI';
% perfIndex = 'ARI';
perfIndex = 'NID';

% Clustering methods
methodNames = {'k-DBA', 'k-Shape', 'k-GAK', 'k-SC', 'k-VARs (oracle)'};

%% Loading and formating results of clusterin precision
resExpr = cell(nPara, 0);

% load and compute results for kVARs
load('../result_expr3/kVARs_result.mat')
for iPara = 1:nPara
    groupKVARs = labelKVARs{iPara};

    precKVARs = [];
    for iExpr = 1:nExpr
        labels = labelExpr{iPara, iExpr};

        groups = groupKVARs(iExpr, :)';
        if ~any(isnan(groups()))
            switch perfIndex
              case 'RI'
                accuracy = perfRI(groups, labels, K);
              case 'NMI'
                accuracy = perfNMI(groups, labels, K);
              case 'ARI'
                accuracy = perfARI(groups, labels, K);
              case 'NID'
                accuracy = 1 - perfNID(groups, labels, K);
              otherwise
                [accuracy, map] = evalAccuracy(groups, labels, K);
            end
        else
            accuracy = nan;
        end
        precKVARs = [precKVARs; accuracy];
    end
    idx = find(contains(methodNames,'k-VARs (oracle)'));
    resExpr{iPara, idx} = precKVARs;
end

% % load and compute results for kVARs (rnd)
% load('result_expr1/kVARs_result_rnd.mat')
% for iPara = 1:nPara
%     groupKVARs = labelKVARs{iPara};

%     precKVARs = [];
%     for iExpr = 1:nExpr
%         labels = labelExpr{iPara, iExpr};

%         groups = groupKVARs(iExpr, :)';
%         if ~any(isnan(groups()))
%             switch perfIndex
%               case 'RI'
%                 accuracy = perfRI(groups, labels, K);
%               case 'NMI'
%                 accuracy = perfNMI(groups, labels, K);
%               case 'ARI'
%                 accuracy = perfARI(groups, labels, K);
%               case 'NID'
%                 accuracy = 1 - perfNID(groups, labels, K);
%               otherwise
%                 [accuracy, map] = evalAccuracy(groups, labels, K);
%             end
%         else
%             accuracy = nan;
%         end
%         precKVARs = [precKVARs; accuracy];
%     end
%     idx = find(contains(methodNames,'k-VARs (rnd)'));
%     resExpr{iPara, idx} = precKVARs;
% end

% load and compute results for kmeansDTW
respath = '../result_expr3/';
for iPara = 1:nPara
    snr = SNRVec(iPara);

    % load raw results
    fname = [respath, 'kMeansDTW_SNR' num2str(snr) '_result.txt'];
    groupExpr = readmatrix(fname, 'FileType', 'text', 'Delimiter', ',');
    groupExpr(~groupExpr) = 2;

    % compute prec
    precKMeansDTW = [];
    for iExpr = 1:nExpr
        labels = labelExpr{iPara, iExpr};
        groups = groupExpr(iExpr, :)';
        if ~any(isnan(groups()))
            switch perfIndex
              case 'RI'
                accuracy = perfRI(groups, labels, K);
              case 'NMI'
                accuracy = perfNMI(groups, labels, K);
              case 'ARI'
                accuracy = perfARI(groups, labels, K);
              case 'NID'
                accuracy = 1 - perfNID(groups, labels, K);
              otherwise
                [accuracy, map] = evalAccuracy(groups, labels, K);
            end
        else
            accuracy = nan;
        end
        precKMeansDTW = [precKMeansDTW; accuracy];
    end
    idx = find(contains(methodNames,'k-DBA'));
    resExpr{iPara, idx} = precKMeansDTW;
end

% load and compute results for kshape
respath = '../result_expr3/';
for iPara = 1:nPara
    snr = SNRVec(iPara);

    % load raw results
    fname = [respath, 'kShape_SNR' num2str(snr) '_result.txt'];
    groupExpr = readmatrix(fname, 'FileType', 'text', 'Delimiter', ',');
    groupExpr(~groupExpr) = 2;

    % compute prec
    precKShape = [];
    for iExpr = 1:nExpr
        labels = labelExpr{iPara, iExpr};
        groups = groupExpr(iExpr, :)';
        if ~any(isnan(groups()))
            switch perfIndex
              case 'RI'
                accuracy = perfRI(groups, labels, K);
              case 'NMI'
                accuracy = perfNMI(groups, labels, K);
              case 'ARI'
                accuracy = perfARI(groups, labels, K);
              case 'NID'
                accuracy = 1 - perfNID(groups, labels, K);
              otherwise
                [accuracy, map] = evalAccuracy(groups, labels, K);
            end
        else
            accuracy = nan;
        end
        precKShape = [precKShape; accuracy];
    end
    idx = find(contains(methodNames,'k-Shape'));
    resExpr{iPara, idx} = precKShape;
end

% load and compute results for kernel kmeans
respath = '../result_expr3/';
for iPara = 1:nPara
    snr = SNRVec(iPara);

    % load raw results
    fname = [respath, 'kMeansGAK_SNR' num2str(snr) '_result.txt'];
    groupExpr = readmatrix(fname, 'FileType', 'text', 'Delimiter', ',');
    groupExpr(~groupExpr) = 2;

    % compute prec
    precKMeansGAK = [];
    for iExpr = 1:nExpr
        labels = labelExpr{iPara, iExpr};
        groups = groupExpr(iExpr, :)';
        if ~any(isnan(groups()))
            switch perfIndex
              case 'RI'
                accuracy = perfRI(groups, labels, K);
              case 'NMI'
                accuracy = perfNMI(groups, labels, K);
              case 'ARI'
                accuracy = perfARI(groups, labels, K);
              case 'NID'
                accuracy = 1 - perfNID(groups, labels, K);
              otherwise
                [accuracy, map] = evalAccuracy(groups, labels, K);
            end
        else
            accuracy = nan;
        end
        precKMeansGAK = [precKMeansGAK; accuracy];
    end
    idx = find(contains(methodNames,'k-GAK'));
    resExpr{iPara, idx} = precKMeansGAK;
end

% load and compute results for k-SC
load('../result_expr3/kSC_result.mat')
for iPara = 1:nPara
    groupKSC = labelKSC{iPara};

    precKSC = [];
    for iExpr = 1:nExpr
        labels = labelExpr{iPara, iExpr};

        groups = groupKSC(iExpr, :)';
        if ~any(isnan(groups()))
            switch perfIndex
              case 'RI'
                accuracy = perfRI(groups, labels, K);
              case 'NMI'
                accuracy = perfNMI(groups, labels, K);
              case 'ARI'
                accuracy = perfARI(groups, labels, K);
              case 'NID'
                accuracy = 1 - perfNID(groups, labels, K);
              otherwise
                [accuracy, map] = evalAccuracy(groups, labels, K);
            end
        else
            accuracy = nan;
        end
        precKSC = [precKSC; accuracy];
    end
    idx = find(contains(methodNames,'k-SC'));
    resExpr{iPara, idx} = precKSC;
end


%% Group boxplots

% xlabel for figures
nMethod = length(methodNames);
xAxisNames = strread(num2str(SNRVec), '%s')';
xLabelName = 'SNR (dB)';
switch perfIndex
  case 'RI'
    yLabelName = 'RI';
  case 'NMI'
    yLabelName = 'NMI';
  case 'ARI'
    yLabelName = 'ARI';
  case 'NID'
    yLabelName = '1 - NID (aka. NMI_{max})';
  otherwise
    yLabelName = 'clustering precision';
end

% group data
prec = [];
for k = 1:nPara
    for j = 1:nMethod
        prec= [prec resExpr{k,j}];
    end
end

% positioning variables for boxplot
innerSpace = .5;
pairSpace = innerSpace*(nMethod-1) + 1;
posSingle = (1:pairSpace:1+pairSpace*(length(xAxisNames)-1))';
% positions = [positions positions+innerSpace positions+2*innerSpace positions+3*innerSpace positions+4*innerSpace];
positions = [];
for i = 1:nMethod
    positions = [positions posSingle+(i-1)*innerSpace];
end
xTickPositions = reshape(positions', 1, []);
groupPositions = mean(positions, 2);

% setup colormap
colorTable = colormap('lines');
colorSelected = colorTable(1:nMethod, :);
color = repmat(colorSelected, length(xAxisNames), 1);
close;

% boxplotting
fig_hl = figure;
set(fig_hl, 'units', 'inches', 'position', [9.4306 11.2778 8.5 2.2917])

bp_hl = boxplot(prec, 'position', xTickPositions);

xlabel(xLabelName);
ylabel(yLabelName);
set(gca,'xtick', groupPositions);
set(gca,'xticklabel', xAxisNames);

box_hl = findobj(gca,'Tag','Box');
for j=1:length(box_hl)
   patch(get(box_hl(j),'XData'),get(box_hl(j),'YData'), color(j,:),...
         'FaceAlpha',.8);
   set(bp_hl(6,j), 'color', [0 0 0], 'linewidth', 1.5);
end
boxch_hl = get(gca, 'Children');

yl = ylim;
ylim([yl(1) 1.02]);

% draw ref. line to separate groups of boxes
hold on
xl = xlim;
numsep = length(SNRVec);
for k = 1:numsep-1
    plot([1 1]*(xl(2)/numsep*k+innerSpace), ylim, 'k:')
end
hold off

% % set legend
% lg_hl = legend(boxch_hl(1:nMethod), methodNames);
% set(lg_hl, 'location', 'eastoutside');

%% export in pdf
pos = get(fig_hl,'Position');
set(fig_hl,'PaperPositionMode','Auto','PaperUnits','Inches', ...
           'PaperSize',[pos(3), pos(4)])
figname = ['TSP_boxplot_vSNR_' perfIndex '.pdf'];
print(fig_hl, figname, '-dpdf', '-r0')
