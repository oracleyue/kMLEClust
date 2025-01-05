function valRI = perfRI(labels_res, labels_gt, K)
% PERFRI computes the Rand Index that evaluates the clustering
% performance.
%
% INPUT:
%   labels_res   :   the results of clustering labels
%                    labels_res(i) denotes the cluster label of the i-th sample
%   labels_gt    :   the ground truth of clustering labels
%                    labels_gt(i) denotes the cluster label of the i-th sample
%   K            :   number of clusters
%
% OUTPUT:
%   valRI     :   scalar in [0, 1]; Rand Index, 1 (best)

% Copyright (c) 2020, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last modified on 28 May 2021


% Dimension checking
assert(length(labels_gt) == length(labels_res), ...
       'The lengths of the ground truth and the result must be the same.');
N = length(labels_gt);

% Reshape results
Yres = zeros(K, N);
Ygt = zeros(K, N);
for n = 1:N
    Yres(labels_res(n), n) = 1;
    Ygt(labels_gt(n), n) = 1;
end

% Computing RandIndex
[mYs, nYs] = size(Yres);
[mYt, nYt] = size(Ygt);
TP=0; TN=0; FN=0; FP=0;
for i = 1:nYs
    for j = i+1:nYs
        n1 = sum(Yres(:,i) == Yres(:,j));
        n2 = sum(Ygt(:,i) == Ygt(:,j));
        if ((n1 == mYs) && (n2 == mYt))
            TP = TP + 1;
        elseif((n1 < mYs) && (n2 < mYt))
            TN = TN + 1;
        elseif((n1 == mYs) && (n2 < mYt))
            FN = FN + 1;
        else
            FP = FP + 1;
        end
    end
end

% valRI = (TP + TN) / (TP + TN + FN + FP);
valRI = (TP + TN) / (N*(N-1)/2);