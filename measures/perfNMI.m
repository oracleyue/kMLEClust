function valNMI = perfNMI(labels_res, labels_gt, K)
% PERFNMI computes the normalized mutual information (NMI) that
% evaluates the clustering performance.
%
% INPUT:
%   labels_res   :   the results of clustering labels
%                    labels_res(i) denotes the cluster label of the i-th sample
%   labels_gt    :   the ground truth of clustering labels
%                    labels_gt(i) denotes the cluster label of the i-th sample
%   K            :   number of clusters
%
% OUTPUT:
%   valNMI       :   scalar in [0, 1]; NMI, 1 (best)

% Copyright (c) 2020, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last modified on 28 May 2021


% dimension checking
assert(length(labels_gt) == length(labels_res), ...
       'The lengths of the ground truth and the results must be the same.');
N = length(labels_gt);

% Reshape results
Yres = zeros(K, N);
Ygt = zeros(K, N);
for n = 1:N
    Yres(labels_res(n), n) = 1;
    Ygt(labels_gt(n), n) = 1;
end

% compute numerator of NMI
num = 0;
for i = 1:K
    for j = 1:K
        Nij = sum(Ygt(i,:) & Yres(j,:));
        NGi = sum(Ygt(i,:));
        NAj = sum(Yres(j,:));
        if ~Nij
            ; % num += 0
        else
            num = num + Nij * log((N*Nij) / (NGi * NAj));
        end
    end
end

% compute denominator of NMI
sumNG = 0;
sumNA = 0;
for i = 1:K
    NGi = sum(Ygt(i,:));
    NAi = sum(Yres(i,:));
    if ~NGi
        ; % +=0
    else
        sumNG = sumNG + NGi * log(NGi / N);
    end
    if ~NAi
        ; % +=0
    else
        sumNA = sumNA + NAi * log(NAi / N);
    end
end
den = sqrt(sumNG * sumNA);      % USSL and other clustering papers all use this index
% den = - (sumNG + sumNA) / 2;  % by def. of NMI, it seems this one

% return NMI index
if den == 0
    valNMI = 0;
else
    valNMI = num / den;
end
