function [labels, centers] = kmeansVARs(tsdata, K, p, varargin)
% KMEANSVARS applies the k-means on VAR models of each time series
% for time-series clustering. It is two-step approach, the first to
% estimate VAR of each time series, then applying k-means on the
% data matrix constructed by all VAR model parameters.
%
% INPUT:
%   tsdata    :   N-dim cells, each element is an (mxT) matrix
%                 N - number of samples, each of which is a multivariate
%                     time series
%                 m - dimention of each point in time series
%                 T - length of time series
%   K         :   positive integer
%                 number of groups/clusters
%   p         :   positive integer
%                 model/lag order of VARs, p1 = ... = pK = p
%                 (not supported) K-dim vector of positive integers
%                 (not supported) model order of VARs, [p1, ..., pK]
%
% OUTPUT:
%   labels   :   (N x 1) vector of integers in {1, ..., K}
%                the i-th value tells the group label of the i-th sample
%   centers  :   k-means centers for AR parameters
%
% Copyright (c) 2020, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last update on 06 Dec 2023


% Argument Parsing
% validation functios
validInteger = @(x) isnumeric(x) && isscalar(x) && ~rem(x,1);
validP = @(x) isvector(x) && ~any(rem(x,1)) && (x > 0);
% setup parser
parser = inputParser;
parser.addRequired('tsdata', @iscell);
parser.addRequired('K', validInteger);
parser.addRequired('p', validP);
% parsing
parser.parse(tsdata, K, p, varargin{:});

% dimension checking
N = length(tsdata);
[m, T] = size(tsdata{1});
assert((T-p) > (1+m*max(p)), ...
       'The length of time series is not sufficient.');

% ================================================
% Algorithm: ts-specific VAR Estimation + k-Means
% ================================================

% Precompution
mdlVARs = cell(2, N);  % 1: Theta, 2: Omega
Thetas = cell(1, N);
Omegas = cell(1, N);
parfor n = 1:N
    % % navie VAR estimation
    % [Xnk, Ynk] = varDataSetup(tsdata{n}, p);
    % [Qnk, Rnk] = qr(Xnk, 0);
    % Theta = Rnk \ (Qnk' * Ynk);
    % Enk = Ynk - Xnk * Theta;
    % [Unk, Vnk] = qr(Enk, 0);
    % Omega = 1/(T-p) * Vnk'*Vnk;

    % use matlab built-in
    mdl = varm(m, p);
    estMdl = estimate(mdl, tsdata{n}');
    Theta = [estMdl.Constant cell2mat(estMdl.AR)];
    Omega = estMdl.Covariance;

    % saving
    Thetas{n} = Theta;
    Omegas{n} = Omega;
end
mdlVARs(1, :) = Thetas;
mdlVARs(2, :) = Omegas;

% perform k-means on VAR paremeters
param = [];
for n = 1:N
%     param = [param; reshape(mdlVARs{1, n}, 1, [])];
    % param = [param; reshape(mdlVARs{2, n}, 1, [])];
    param = [param; reshape(mdlVARs{1, n}, 1, []), reshape(mdlVARs{2, n}, 1, [])];
end
[labels, centers] = kmeans(param, K);

end % END of Interface


% ================================================================
% Local Functions
% ================================================================

function [X, Y] = varDataSetup(data, p)
% Set up the data matrices X, Y for VAR estimation.
% Notes:
%  - data: (mxT) matrix, row indexing variable dim, and column
% indexing time.
%  - p: integer, the model order of VAR.
%  - X: (T-p)x(1+mp) matrix, \mathbf{X}_{nk}
%  - Y: (T-p)xmp matrix, \mathbf{Y}_{nk}

    [m, T] = size(data);
    X = []; Y = [];

    for t = p+1:T
        Yt = data(:, t);
        Xkt = [1; reshape(fliplr(data(:, (t-p):(t-1))), [], 1)];

        % intermediate, pending to be transposed
        X = [X, Xkt];
        Y = [Y, Yt];
    end
    X = X';
    Y = Y';
end % END of varDataSetup
