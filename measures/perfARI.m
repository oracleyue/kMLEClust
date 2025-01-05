function valRI = perfARI(res, gt, K)
% PERFRI computes the adjusted Rand Index (ARI) that evaluates the
% clustering performance.
%
% INPUT:
%   res       :   the results of clustering labels
%                 res(i) denotes the cluster label of the i-th sample
%   gt        :   the ground truth of clustering labels
%                 gt(i) denotes the cluster label of the i-th sample
%   K         :   number of clusters
%
% OUTPUT:
%   valRI     :   scalar in [0, 1]; Rand Index, 1 (best)
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

% Compute components
N00 = 0; N11 = 0; N01 = 0; N10 = 0;
for i = 1:N
    for j = 1:N
        % data pair (i,j): res for U, gt for V in paper
        % Reference: Vinh, N. X., Epps, J. and Bailey, J. (2010) ‘Information theoretic measures for clusterings comparison: Variants, properties, normalization and correction for chance’, The Journal of Machine Learning Research. JMLR. org, 11, pp. 2837–2854.
        if (res(i) == res(j)) & (gt(i) == gt(j))
            N11 = N11 + 1;
        elseif  (res(i) ~= res(j)) &  (gt(i) ~= gt(j))
            N00 = N00 + 1;
        elseif (res(i) == res(j))  &  (gt(i) ~= gt(j))
            N01 = N01 + 1;
        elseif  (res(i) ~= res(j)) & (gt(i) == gt(j))
            N10 = N10 + 1;
        else
            error('Heartbleed bug, debug!')
        end
    end
end

% Evaluate ARI
num = 2 * (N00 * N11 - N01 * N10);
den = (N00 + N01) * (N01 + N11) + (N00 + N10) * (N10 + N11);
valRI = num / den;
