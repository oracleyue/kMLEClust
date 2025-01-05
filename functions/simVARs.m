function [data, mdls, labels] = simVARs(m, p, T, N, K, spec)
% SIMVARS randomly generates the specified number of VAR models and
% simulates each model the requested amount of times.
%
% INPUT:
%   m   :   positive integer
%           Dimension of variates, i.e., output Yt is m-dimensional
%   p   :   positive integer
%           Model order (or lag order) of VAR models, i.e. VAR(p)
%           OR, K-dim vector of positive integers
%           Model orders specifying K number of VAR models, i.e. VAR(p1), ... VAR(pK)
%   T   :   positive integer
%           Length of time serie to be generated
%   N   :   positive integer (default: 1)
%           if spec is 'total' or unspecified, N denotes the total number of
%              time series to be generated;
%           if spec is 'per_cluster', N denotes the number of time series per cluster.
%   K   :   positive integer (default: 1)
%           Number of VAR(p) models to be generated
%   spec:   string, 'total' or 'per_cluster'
%           Control the meaning of argument "N"
%
% OUTPUT:
%   data      :   N-dim cell, OR (mxT) matrix if N = 1
%                 Each element is an (mxT) matrix that represents a
%                 multivariate time series of dim m and time length T
%   mdls      :   K-dim cell, OR "varm" model if K = 1
%                 Each element is a "varm" model that represents m-dim VAR(p)
%   labels    :   N-dim vector consists of integers in [1:K]
%                 It specifies a time series is simulated by which VAR model,
%                 e.g., labels(5) = 2 implies the 5th time series is
%                 generatated by the 2nd VAR model in "mdls".
%
% Examples:
%  1. To randomly generate a 3-dim VAR(2) and simulate a time series
%     of length 100:
%     [data, mdl] = simVARs(3, 2, 100);
%
%  2. To randomly generate two 3-dim VAR(2) models, and
%     simulate 20 number of time series of length 100 for clustering
%     [data, mdl, labels] = simVARs(3, 2, 100, 2, 20);

% Copyright (c) 2020, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last update on 28 Aug 2020


% Prerequisites
% - Econometrics Toolbox
%   Functions: varm, simulate

if ~license('test','econometrics_toolbox')
    error('This function requires Econometrics Toolbox, incl. varm, simulate.')
end


% Debug flags
plotFlag = 0;


% Parsing arguments
if nargin < 3
    error('The arguments "m,p,T" are required.')
elseif nargin == 3
    K = 1; N = 1; spec = 'total';
elseif nargin == 4
    K = 1; spec = 'total';
elseif nargin == 5
    spec = 'total';
elseif nargin == 6
    switch spec
      case 'total'
        ;
      case 'per_cluster'
        Nc = N;
        N = Nc * K;
    end
end


% Random generation of stable VAR(pk) models
if K > 1
    mdls = cell(K, 1);
end

% model orders (p1,...,pK) for K VAR models
if isscalar(p)
    p = ones(K, 1)*p;
end

% set up VAR(pk) models
for k = 1:K
    % Theta0, or varm.constant
    const = randn(m, 1);

    % Thetas, or varm.AR
    lag = p(k);
    rt = rand(m, lag) + 2;    % rt >> 1 vibrates more sharply
    lambdaM = zeros(m, lag);
    for i = 1:m
        pl = poly(rt(i, :));
        pl = fliplr(pl);
        pl = pl ./ pl(1);
        lambdaM(i, :) = -pl(2:end);
    end
    Thetas = [];
    temp = randn(m, m);
    [LTheta, ~] = qr(temp);
    for l = 1:p
        eigTheta = lambdaM(:,l);
        Theta = LTheta' * diag(eigTheta) * LTheta;
        Thetas = [Thetas, {Theta}];
    end

    % Oemga, or var.covariance
    LOmega = randn(m, m); % it is singular with prob. 0.
                          % if like to ensure non-singularity, use the following.
    % while 1
    %     LOmega = randn(m, m);
    %     LOmega = tril(LOmega);
    %     if any(~diag(LOmega))
    %         continue
    %     end
    %     break
    % end
    Omega = LOmega' * LOmega;
    % being nomalized if needed
    % Omega = Omega / norm(Omega, 'fro');

    % set up varm object
    mdl = varm('constant', const, 'AR', Thetas, 'covariance', Omega);

    % save for returns
    if K > 1
        mdls{k} = mdl;
    else
        mdls = mdl;
    end
end


% Simulate VAR(p1),...,VAR(pK) models
if N > 1
    data = cell(N, 1);
end

% labels that specify which VAR model is used to generate the n-th signal
switch spec
  case 'total'
    labels = randi(K, [1, N])';
  case 'per_cluster'
    labels = kron((1:K)', ones(Nc,1));
    labels = labels(randperm(N));
end

% simulation
for n = 1:N
    idxVAR = labels(n);
    if K > 1
        mdl = mdls{idxVAR};
    else
        mdl = mdls;
    end
    ts = simulate(mdl, T);

    % save for returns
    if N > 1
        data{n} = ts';
    else
        data = ts';
    end
end


% Plot for debugging or checking
if plotFlag
    idxN = randi(N);
    if N > 1
        ts = data{idxN};
    else
        ts = data;
    end
    if m > 4
        idxChan = sort(randi(m, [1 4]));
    else
        idxChan = 1:m;
    end

    fighl = figure;
    legstr = [];
    hold on
    if m > 4
        for i = 1:4
            plot(1:T, ts(idxN(i),:), 'o-');
            legstr = [legstr, {sprintf('variate #%d', idxN(i))}];
        end
    else
        for i = 1:m
            plot(1:T, ts(i,:), 'o-');
            legstr = [legstr, {sprintf('variate #%d', i)}];
        end
    end
    hold off
    grid on
    xlabel('time');
    legend(legstr);
    title(sprintf('Sample plot of time series #%d from VAR #%d', ...
                  idxN, labels(idxN)));
end