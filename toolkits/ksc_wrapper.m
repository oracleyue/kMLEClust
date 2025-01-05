function labels = ksc_wrapper(data, K)
% KSC_WRAPPER is a wrapper function uses the k-SC algorithm from
% "k-SC" toolbox to clustering multivariate time series.
%
% INPUT:
%   data      :   N-dim cell
%                 Each element is an (mxT) matrix that represents a
%                 multivariate time series of dim m and time length T
%   K         :   number of cluster
% OUTPUT:
%   labels   :   (N x 1) vector of integers in {1, ..., K}
%                the i-th value tells the group label of the i-th sample

% Copyright (c) 2021, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last update on 28 May 2021

% Idea: since the ksc_toy() can only cluster univariate time series, we
% use treate each channel of a vector time series as a separate one; and
% the resultant cluster of a vector time series is determined by the
% most occurance among its channels.

% dimensions
N = length(data);
[m, T] = size(data{1});
labels = zeros(N, 1);

% concatenate all time series
Udata = [];
for n = 1:N
    Udata = [Udata; data{n}];
end

% clustering by ksc
Ulabels = ksc_toy(Udata, K);

% determine cluster of each vector ts
for n = 1:N
    startIdx = (n-1)*m+1;
    endIdx = n*m;
    labels(n) = mode(Ulabels(startIdx:1:endIdx));
end
