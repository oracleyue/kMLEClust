function valNID = perfNID(res, gt, K)
% PERFRI computes the Normalized Information Distance (NID) that
% evaluates the clustering performance.
%
% INPUT:
%   res       :   the results of clustering labels
%                 res(i) denotes the cluster label of the i-th sample
%   gt        :   the ground truth of clustering labels
%                 gt(i) denotes the cluster label of the i-th sample
%   K         :   number of clusters
%
% OUTPUT:
%   valNID     :   scalar in [0, 1]; Rand Index, 1 (best)
%
% REFERENCE:
% Vinh, N. X., Epps, J. and Bailey, J. (2010) ‘Information theoretic
% measures for clusterings comparison: Variants, properties,
% normalization and correction for chance’, The Journal of Machine
% Learning Research. JMLR. org, 11, pp. 2837–2854.

% Copyright (c) 2020, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last modified on 26 Mar 2022


% Dimension checking
assert(length(gt) == length(res), ...
       'The lengths of the ground truth and the result must be the same.');
N = length(gt);

% Reformating U (res) and V (gt)
% U(i,:) denotes the samples in cluster U_i
U = zeros(K, N); V = zeros(K, N);
for n = 1:N
    U(res(n), n) = true;
    V(gt(n), n) = true;
end

% Compute components n_ij, a_i, b_j in paper
Ntbl = zeros(K, K);
aList = zeros(K, 1);
bList = zeros(K, 1);
for i = 1:K
    for j = 1:K
        % cluster i in res; cluster j in gt
        Ntbl(i,j) = sum(and(U(i,:), V(j,:)));
    end
    aList(i) = sum(U(i,:));
    bList(i) = sum(V(i,:));
end

% Compute mutual information indexes
HU = 0; HV = 0; IUV = 0;
for i = 1:K
    for j = 1:K
        num = Ntbl(i,j) / N;
        den = aList(i) * bList(j) / N^2;
        if num == 0
            IUV = IUV + 0;  % Informatics: 0log(0) = 0
        else
            IUV = IUV + num * log(num / den);
        end
    end
    if aList(i) ~= 0  % 0log0 = 0, then HU no change
        HU = HU - (aList(i)/N) * log(aList(i)/N);
    end
    if bList(i) ~= 0  % 0log0 = 0, then HU no change
        HV = HV - (bList(i)/N) * log(bList(i)/N);
    end
end

% Evaluate NID
valNID = 1 - IUV / max([HU, HV]);
