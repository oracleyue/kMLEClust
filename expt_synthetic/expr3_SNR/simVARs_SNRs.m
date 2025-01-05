function [dataList, mdlList, labelList, SNRrList] = simVARs_SNRs(m, p, T, SNRList, N, K, spec)
% SIMVARS_SNRs randomly generates the specified number of VAR models and
% simulates each model the requested amount of times, with different levels of noises.
% It returns an array of models that corresponds to the list of SNRs.
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
%   SNRList :   vector of integers
%           Signal-noise-ratio, 10*log(P_signal/P_noise), unit: dB
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
%   dataList   :   cell array of "data"
%   mdlList    :   cell array of "mdls"
%   labelList  :   cell array of "labels"
%   SNRrTbl       :   actural SNR values for all K models
%   data      :   N-dim cell, OR (mxT) matrix if N = 1
%                 Each element is an (mxT) matrix that represents a
%                 multivariate time series of dim m and time length T
%   mdls      :   K-dim cell, OR "varm" model if K = 1
%                 Each element is a "varm" model that represents m-dim VAR(p)
%   labels    :   N-dim vector consists of integers in [1:K]
%                 It specifies a time series is simulated by which VAR model,
%                 e.g., labels(5) = 2 implies the 5th time series is
%                 generatated by the 2nd VAR model in "mdls".

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
    SNRList = 10;
    K = 1; N = 1; spec = 'total';
elseif nargin == 4
    K = 1; N = 1; spec = 'total';
elseif nargin == 5
    K = 1; spec = 'total';
elseif nargin == 6
    spec = 'total';
elseif nargin == 7
    switch spec
      case 'total'
        ;
      case 'per_cluster'
        Nc = N;
        N = Nc * K;
    end
end

% Keep all data corresponding to SNRList
nSNR = length(SNRList);
if nSNR > 1
    dataList = cell(nSNR, 1);
    mdlList = cell(nSNR, 1);
    labelList = cell(nSNR, 1);
    SNRrTbl = cell(K, 1);   % save SNR of every VAR model
end

% Random generation of stable VAR(pk) models
mdls = cell(K, nSNR);

% model orders (p1,...,pK) for K VAR models
if isscalar(p)
    p = ones(K, 1)*p;
end

% set up VAR(pk) models
for k = 1:K
    % Oemga, or var.covariance
    LOmega = randn(m, m); % it is singular with prob. 0.
    Omega = LOmega' * LOmega;
    % being nomalized if needed
    Omega = Omega / norm(Omega, 'fro');
    % for computing SNR
    maxEigenOm = abs(eigs(Omega,1));

    % Theta0, or varm.constant
    const = randn(m, 1);

    % Roots for building Thetas (varm.AR)
    lag = p(k);
    rt = rand(m, lag);
    for i = 1:m
        rt(i,randi(lag)) = 10^(-randi(4));  % rt close to 1 for larger SNR
    end
    rt = rt + 1;   % ensure rts > 1

    % Building table for alpha vs. SNR for table search for SNRList
    alphaTable = 1:0.001:5;
    for iAlph = 1:length(alphaTable)

        % scaling roots
        alpha = alphaTable(iAlph);
        roots = alpha * rt;

        % building As from roots
        As = buildVARPara(roots);

        % SNR computation: snr = 10 log(rho(Gamma0,Omega))
        % rho(X,M): Reimainian distance between p.d. X and M
        Gamma0 = cov0VAR(As, Omega);
        SNRTable(iAlph) = 10 * log10(ReimainDist(Gamma0, Omega));
    end

    disp('Alpha-SNR table building: done.')

    for iSNR = 1:nSNR
        [val, idx] = min(abs(SNRTable - SNRList(iSNR)));
        alphaList(iSNR) = alphaTable(idx);
        SNRrList(iSNR) = SNRTable(idx);
    end

    % plotting for debug
    if plotFlag
        figure;
        plot(alphaTable, SNRTable, '-*')
        xlabel('scaling factor')
        ylabel('SNR')
        title(sprintf('For SNR of cluster K=%d', k))
    end
    % disp actual SNRs for models
    fprintf('SNR values for model-%d:',k);
    disp(SNRrList)

    % set up and save varm objects
    for iSNR = 1:nSNR
        alpha = alphaList(iSNR);
        As = buildVARPara(alpha*rt);
        mdl = varm('constant', const, 'AR', As, 'covariance', Omega);
        mdls{k,iSNR} = mdl;
    end
    SNRrTbl{k} = SNRrList;
end

% Simulating for nSNR number of model sets
for iSNR = 1:nSNR
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
            mdl = mdls{idxVAR, iSNR};
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

    % saving for this SNR
    mdlList{iSNR} = mdls(:, iSNR);
    dataList{iSNR} = data;
    labelList{iSNR} = labels;

    % plot for debugging or checking
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
end

end % END of all

% ================================================================
% Local Functions
% ================================================================
function As = buildVARPara(rt)
    % Contructing VAR parameters A(1),...,A(p) from roots.

    [m, p] = size(rt);

    lambdaM = zeros(m, p);
    for i = 1:m
        pl = poly(rt(i, :));
        pl = fliplr(pl);
        pl = pl ./ pl(1);
        lambdaM(i, :) = -pl(2:end);
    end
    As = [];
    temp = randn(m, m);
    [LTheta, ~] = qr(temp);
    for l = 1:p
        eigTheta = lambdaM(:,l);
        Theta = LTheta' * diag(eigTheta) * LTheta;
        As = [As, {Theta}];
    end
end

function Gamma0 = cov0VAR(As, Omega)
    % Compute the autocovariance of VAR at 0 via Yuler-Walker equation.

    % dimensions
    p = length(As);
    m = size(As{1}, 1);

    % construct Yt = nu + A Yt-1 + Ut  (Lutkepohl, 2005; pp.15)
    Arow = cat(2, As{1:end});
    A = [Arow; [eye(m*(p-1)), zeros(m*(p-1),m)]];
    OmegaU = blkdiag(Omega, zeros(m*(p-1)));

    % solving for Gamma: Gamma = A*Gamma*A' + OmegaU
    % 1. naive inverse method:
    % Akron = eye(m*p*m*p) - kron(A, A);
    % GammaVec = Akron \ reshape(OmegaU, [], 1);
    % Gamma = reshape(GammaVec, m*p, m*p);
    % 2. Lyapunov eq.
    Gamma = dlyap(A, OmegaU);

    % extract autocov at 0
    Gamma0 = Gamma(1:m, 1:m);
end

% function maxval = maxEigVAR(As)
%     % Compute the maximum eigenvalue of VAR parameter matrix.

%     % dimensions
%     p = length(As);
%     m = size(As{1}, 1);

%     % construct Yt = nu + A Yt-1 + Ut  (Lutkepohl, 2005; pp.15)
%     Arow = cat(2, As{1:end});
%     A = [Arow; [eye(m*(p-1)), zeros(m*(p-1),m)]];

%     maxval = max(eig(A));
% end

function dist = ReimainDist(Gamma, Omega)
    % Compute the Reimainian distance between Gamma and Omega.

    % Way 1: using matrix square root
    % OmRoot = sqrtm(Omega);
    % eigVal = eig(OmRoot \ Gamma / OmRoot);

    % Way 2: using eigen-decomposition
    [V, D] = eig(Omega);
    Drt = sqrt(inv(D));    % fast since D is diagonal
    eigVal = eig(Drt * (V \ Gamma * V) * Drt);

    % Reimainian distance: rho(X,M) = sqrt(sum(ln^{2}(M^{-1/2} X M^{-1/2})))
    dist = sqrt(sum(log(eigVal).^2));
end
