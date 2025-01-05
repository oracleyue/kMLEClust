% load the data file first.
% load('Wafer.mat')

% read data info
train = mts.train;
test = mts.test;
trainlabels = mts.trainlabels;
testlabels = mts.testlabels;

numTrain = length(train);
numTest = length(test);
numCluster = max(trainlabels);

dimVar = size(train{1}, 1);
dimTsLenTrain = [];
dimTsLenTest = [];

for k = 1:length(train)
    dimTsLenTrain = [dimTsLenTrain size(train{k}, 2)];
end
for k = 1:length(test)
    dimTsLenTest = [dimTsLenTest size(test{k}, 2)];
end

% summary
fprintf('#Samples in the training set: %d\n', numTrain);
fprintf('#Samples in the test set: %d\n', numTest);
fprintf('#Clusters: %d\n', numCluster);
fprintf('Variable dimension: %d\n', dimVar);
fprintf('Lengths of time series (training set): %d-%d\n',...
        min(dimTsLenTrain), max(dimTsLenTrain));
fprintf('Lengths of time series (test set): %d-%d\n',...
        min(dimTsLenTest), max(dimTsLenTest));
