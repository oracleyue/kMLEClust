function [labels, mdl, lld] = kVARs(tsdata, K, p, varargin)
% KVARS applies the k-VARs approach to clustering time-series data,
% which is an extended version of k-means for time series.
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
%   [argument-pairs]
%   "tol"     :   positive float
%                 tolerance, default 1e-8
%   "maxIter" :   positive integer
%                 maximum number of EMiterations, default 100
%   "initType":   char: 'random' (default), 'kmeans'
%                 NOTE: if "init" provided, this argument is obsoleted.
%                 'random' - randomly choose signals to id. VARs for init
%                 'kmeans', 'kmeansC' - clustering N VAR paras for init
%                 'user' - user-specified by "init" option
%   "init"    :   if K-dim vector of different integers ranging in 1:N
%                    user-specified signals for initializing K VAR models
%                    used only for "initType = user" (automatically set "initType")
%                 if an integer, number of tirals for random initialization
%                    used only for "initType = random" (automatically set "initType")
%
% OUTPUT:
%   labels   :   (N x 1) vector of integers in {1, ..., K}
%                the i-th value tells the group label of the i-th sample
%   mdl      :   struct
%                mdl.Alpha: K-dim vector of doubles in (0,1)
%                mdl.Theta: K-dim cell of matrices as VAR parameters,
%                           each element represents [Theta0, Theta1, ..., Theta_p]
%                mdl.Omega: K-dim cell of matrices as VAR covariances
%                mdl.Tau: (NxK) matrix of {0,1}, each row is 1-hot
%                mdl.Phi: (NxK) matrix, prediction error of k-th VAR for n-th time series
%   lld      :   double; the value of log-likelihood (up to a certain constant)
%
% Examples:
%   labels = vecKARs(tsdata, 2, 5);
%   labels = vecKARs(tsdata, 2, [2 2 3 5 5]);
%   labels = vecKARs(tsdata, 2, 5, 'tol', 1e-6, 'maxIter', 200);
%   labels = vecKARs(tsdata, 2, 5, 'init', [1 4]);

% Copyright (c) 2020, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last update on 25 Aug 2020


% Argument Parsing

% default values
defaultTol = 1e-8;
defaultMaxIter = 100;
defaultInit = [];
defaultInitType = 'random';
% defaultStopType = 'zeroParaChanges';
% validation functios
validInteger = @(x) isnumeric(x) && isscalar(x) && ~rem(x,1);
validP = @(x) isvector(x) && ~any(rem(x,1)) && (x > 0);
validTol = @(x) isnumeric(x) && isscalar(x) && (x < 1) && (x > 0);
validMaxIter = @(x) validInteger(x);
validInit = @(x) isnumeric(x) && isvector(x);
validInitType = @(x) any(strcmp(x, {'random', 'kmeans', 'user'}));
% validStopType = @(x) any(strcmp(x, {'zeroCostChanges', 'zeroParaChanges'}));
% setup parser
parser = inputParser;
parser.addRequired('tsdata', @iscell);
parser.addRequired('K', validInteger);
parser.addRequired('p', validP);
parser.addParameter('tol', defaultTol, validTol);
parser.addParameter('maxIter', defaultMaxIter, validMaxIter);
parser.addParameter('init', defaultInit, validInit);
parser.addParameter('initType', defaultInitType, validInitType);
% parser.addParameter('stopType', defaultStopType, validStopType);
% parsing
parser.parse(tsdata, K, p, varargin{:});
tol = parser.Results.tol;
maxIter = parser.Results.maxIter;
initType = parser.Results.initType;
% stopType = parser.Results.stopType;
init = sort(parser.Results.init);
initArgCheck = (length(init) == 1 && any(strcmp(initType, 'random'))) || ...
               (length(init) > 1  && any(strcmp(initType, 'user'))) || ...
               (length(init) == 0);
assert(initArgCheck, 'The compatibility of arguments "init" and "initType" fails!');
if length(init) == K
    initType = 'user';
    idxInit = init;
elseif length(init) == 1
    initType = 'random';
    numInit = init;
    idxInit = [];
elseif isempty(init)
    idxInit = [];
else
    error('Argument "init" is not valid!')
end
if isempty(init) && strcmp(initType, 'random')
    numInit = 1;
end

% dimension checking
N = length(tsdata);
[m, T] = size(tsdata{1});
assert((T-p) > (1+m*max(p)), ...
       'The length of time series is not sufficient.');

% check illegal setting of "idxInit"
checkIdxInit = isempty(idxInit) || ...
    (length(idxInit) == K) && (max(idxInit) <= N) && ...
    (length(idxInit) == length(unique(idxInit)));
assert(checkIdxInit, 'Argument "init" must be a K-dim vector of unique integers in [1, N].')

% perform k-VARs according to "initType"
switch initType
    case 'random'
        for k = 1:numInit
            % initialize by random K time series
            idxInit = sort(randi(N, [1, K]));
            % runs k-VARs
            [labelsTrial, mdlTrial, lldTrial] = clustKVARs(tsdata, K, p, idxInit, [tol, maxIter]);
            if (k == 1) || (lld < lldTrial)
                labels = labelsTrial;
                mdl = mdlTrial;
                lld = lldTrial;
            end
        end

    case 'kmeans'
        % initialize by k-means of VAR parameters
        if isempty(idxInit)
            paraARs = [];
            for n = 1:N
                paraARs = [paraARs; reshape(mdlPreCell{1, n}, 1, [])];
            end
            [idxKmeans, centKmeans] = kmeans(paraARs, K);
        end
        idxInit = zeros(1, K);
        for k = 1:K
            idxGrpK = find(idxKmeans == k);
            if ~isempty(idxGrpK)
                idx = randi(length(idxGrpK), 1);
                idxInit(k) = idxGrpK(idx);
            else
                msgID = 'kVARs:BadInitialization';
                msgtext = 'Initialization by pre K-means failed. You may try random initialization, or use the argument "init" to specify an initialization.';
                INITE = MException(msgID,msgtext);
                throw(INITE)
            end
        end
        % running k-VARs
        [labels, mdl, lld] = clustKVARs(tsdata, K, p, idxInit, [tol, maxIter]);

    case 'user'
        % running k-VARs
        [labels, mdl, lld] = clustKVARs(tsdata, K, p, idxInit, [tol, maxIter]);
end

end % END of Interface

% ================================================================
% Main Function
% ================================================================

function [labels, mdl, lld] = clustKVARs(tsdata, K, p, idxInit, iterOpt)
    % This is the main function that runs k-VARs to cluster tsdata
    % according to the specified cluster initialization.

    % Dimensions
    N = length(tsdata);
    [m, T] = size(tsdata{1});

    % Algorithm precisions
    tol = iterOpt(1);
    maxIter = iterOpt(2);

    % Precompution
    mdlPreCell = cell(2, N);  % 1: Theta, 2: Omega
    XnkCell = cell(1, N);
    YnkCell = cell(1, N);
    QnkCell = cell(1, N);
    RnkCell = cell(1, N);

    kThetas = cell(1, N);
    kOmegas = cell(1, N);
    parfor n = 1:N
        [Xnk, Ynk] = varDataSetup(tsdata{n}, p);
        [Qnk, Rnk] = qr(Xnk, 0);
        kThetaPre = Rnk \ (Qnk' * Ynk);
        Enk = Ynk - Xnk * kThetaPre;
        [Unk, Vnk] = qr(Enk, 0);
        kOmegaPre = 1/(T-p) * Vnk'*Vnk;

        QnkCell{n} = Qnk;
        RnkCell{n} = Rnk;
        XnkCell{n} = Xnk;
        YnkCell{n} = Ynk;
        kThetas{n} = kThetaPre;
        kOmegas{n} = kOmegaPre;
    end
    mdlPreCell(1, :) = kThetas;
    mdlPreCell(2, :) = kOmegas;

    % Initialization
    Thetas = mdlPreCell(1, idxInit);
    Omegas = mdlPreCell(2, idxInit);
    Alphas = ones(1, K) * (1/K);

    % Expectation Maximization
    countIter = 0;
    while 1
        % Backup prev. parameter before any updating
        prevThetas = Thetas;
        prevOmegas = Omegas;
        prevAlphas = Alphas;

        % E-step
        PhiNK = zeros(N, K);
        TauNK = logical(zeros(N, K));
        cardIK = zeros(1, K);
        parfor n = 1:N
            Xnk = XnkCell{n};
            Ynk = YnkCell{n};
            Qnk = QnkCell{n};
            Rnk = RnkCell{n};

            for k = 1:K
                PhiNK(n,k) = dissimilarityMeasure(Xnk, Ynk, Thetas{k}, Omegas{k}, 'kMLE');
            end
        end
        [~, labels]= min(PhiNK, [], 2);
        for n = 1:N
            TauNK(n, labels(n)) = 1;
        end

        % M-step
        cardIK = sum(TauNK);
        Alphas = cardIK / N;
        parfor k = 1:K
            idxIk = find(TauNK(:, k));

            % EXCEPTION: no signal clustered to this group
            % This is due to bad initialization. However, it is not true that every bad initialization leads to this exception.
            % if isempty(idxIk)
            %     msgID = 'kVARs:InvalidResponsibilities';
            %     msgtext = 'Bad initialization captured. You may use other seed in rng() to try another random initialization, or use the argument "init" to specify an initialization.';
            %     INITE = MException(msgID,msgtext);
            %     throw(INITE)
            % end

            if isempty(idxIk)
                msg = sprintf('No data is assigned to the %d-th cluster. Model updating is skipped.');
                % warning(msg);
            else
                % update Theta
                XX = 0; XY = 0;  % matrix not scalar
                for n = idxIk'
                    Xnk = XnkCell{n};
                    Ynk = YnkCell{n};
                    % Qnk = QnkCell{n};
                    Rnk = RnkCell{n};

                    XX = XX + Rnk' * Rnk;
                    XY = XY + Xnk' * Ynk;
                    % XY = XY + Rnk' * (Qnk' * Ynk);
                end
                Thetas{k} = XX \ XY;

                % update Omega
                EE = 0;  % matrix not scalar
                for n = idxIk'
                    Xnk = XnkCell{n};
                    Ynk = YnkCell{n};

                    Enk = Ynk - Xnk * Thetas{k};
                    [~, Vnk] = qr(Enk, 0);
                    EE = EE + Vnk' * Vnk;
                end
                Omegas{k} = EE / (T-p) / cardIK(k);
            end
        end

        % Stopping conditions: when zero increments of parameter updating
        incTheta = 0; incOmega = 0;
        for k = 1:K
            incTheta = incTheta + norm(prevThetas{k} - Thetas{k}, 'fro');
            incOmega = incOmega + norm(prevOmegas{k} - Omegas{k}, 'fro');
        end
        incAlpha = norm(prevAlphas - Alphas, 2);

        if max([incTheta, incOmega, incAlpha]) < tol || countIter > maxIter
            % Vector k-ARs model parameters
            mdl.Alpha = Alphas;
            mdl.Theta = Thetas;
            mdl.Omega = Omegas;
            % labels: updated in E-step
            mdl.Tau = TauNK;
            mdl.Phi = PhiNK;
            % lld
            lld = - sum(sum(TauNK .* PhiNK));

            % fprintf('k-VARs: %d iterations\n', countIter)
            break
        end

        % update iteration counter
        countIter = countIter + 1;
    end

    % Format output mdl in a better structure
    mdl = mvarFormat(mdl);

end % END of vecKARs

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


function dm = dissimilarityMeasure(X, Y, Theta, Omega, dmType)
    % Calculate the dissimilarity measure D_nk(\theta_k), which is phi(n,k) for
    % vector k-ARs.
    % Notes:
    %  - X: (T-p)x(1+mp) matrix, \mathbf{X}_{nk}, stack of Xkt for t = p+1,...,T.
    %  - Y: (T-p)xm matrix, \mathbf{Y}_{nk}, stack of Yt for t = p+1,...,T.
    %  - dm: positive scalar
    %  - E: (T-p)xm matrix, \mathbf{E}_{nk}, stack of error signals for t = p+1,...,T.
    %  - Theta: VAR parameter of the k-th VAR, transposed, (mp+1)xm.
    %  - Omega: VAR covariance, mxm positive definite matrix.

    if nargin < 5
        dmType = 'kARs';
    end

    % dimensions
    [Tp, m] = size(Y);
    [mp1, mTheta] = size(Theta);
    assert(m == mTheta, 'Dimensions of arguments fail to match!')
    p = (mp1 - 1) / m;
    T = Tp + p;

    % error signal
    E = Y - X * Theta;

    % chol for Omega and check p.d. of Omega
    try
        L = chol(Omega, 'lower');
    catch ME
        switch ME.identifier
            case 'MATLAB:posdef'
                % chol() fails due to a non-positive definite matrix
                msgID = 'User:NonPositiveDefinite';
                msgtext = 'dissimilarityMeasure(): chol() fails due to non-positive definitiness';
                userME = MException(msgID,msgtext);
                throw(userME);
            otherwise
                rethrow(ME);
        end
    end

    % compute dissimilarity measure
    dm = 0;
    for t = 1:size(E, 1)
        LEt = L \ E(t,:)';
        dm = dm + LEt' * LEt;
    end

    % add logdet(Omega) if using kMLE
    switch dmType
        case 'kARs'
             ;
        case 'kMLE'
            dm = dm + (T-p) * 2*sum(log(diag(L)));
    end
end % END of dissimilarityMeasure


function val = logdet(X, method)
    % Reliably compute log(det(X)), where X must be positive definte matrix.

    if nargin < 2
        method = 'chol';
    end
    assert(any(strcmpi({'eig', 'chol', 'det'}, method)), ...
           'The argument "method" must be "eig", "chol" or "det".');

    switch method
        case 'chol'
            % MATLAB 2019a has a severe bug on "chol"!
            % It fails on matrices with minimal eigenvalue larger than 5e-4!
            [L, flag] = chol(X);
            if flag % chol on a non-positive definite matrix
                msgID = 'User:NonPositiveDefinite';
                msgtext = 'logdet(): chol() fails due to non-positive definitiness';
                ME = MException(msgID,msgtext);
                throw(ME);
            end
            val = 2*sum(log(diag(L)));

        case 'eig'
            eigvalX = eig(X);
            if any(eigvalX <= 0)
                msgID = 'User:NonPositiveDefinite';
                msgtext = 'logdet(): eig() gives non-positive eigenvalues';
                ME = MException(msgID,msgtext);
                throw(ME);
            end
            val = log(prod(eigvalX));

        case 'det'
            % This method is not reliable when Omega is large, e.g. dim > 400.
            val = log(det(X));
            if isinf(val)  % det gives Inf value
                msgID = 'User:InfValueReturned';
                msgtext = 'logdet(): det() gives Inf values';
                ME = MException(msgID,msgtext);
                throw(ME);
            end
    end
end % END of logdet


function mdl = mvarFormat(rawmdl)
    % Format output mdl in a better structure, like adding description and
    % transpose Thetas.

    K = length(rawmdl.Alpha);
    [pm1, m] = size(rawmdl.Theta{1});
    p = (pm1 - 1)/m;

    mdl.Description = ...
        sprintf('%d-Mixture of AR-Stationary %d-Dimensional VAR(%d) Model', ...
                K, m, p);
    mdl.Alpha = rawmdl.Alpha;
    for k = 1:K
        mdl.Theta{k} = rawmdl.Theta{k}';
    end
    mdl.Omega = rawmdl.Omega;

    mdl.Tau = rawmdl.Tau;
    mdl.Phi = rawmdl.Phi;

end %END of mvarFormat
