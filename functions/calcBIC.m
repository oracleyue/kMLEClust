function score = calcBIC(data, mdl)
% CALCBIC computes the BIC scores to choose model order and the number of clusters.
%
% INPUT:
%   data    :    N-dim cells, each element is an (mxT) matrix
%                N - number of samples, each of which is a multivariate
%                    time series
%                m - dimention of each point in time series
%                T - length of time series
%   mdl     :    struct
%                mdl.Alpha: K-dim vector of doubles in (0,1)
%                mdl.Theta: K-dim cell of matrices as VAR parameters,
%                           each element represents [Theta0, Theta1, ..., Theta_p]
%                mdl.Omega: K-dim cell of matrices as VAR covariances
%
% OUTPUT:
%   score   :    positive real value
%

% Copyright (c) 2021, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last update on 12 Feb 2021


% Dimension parsing
N = length(data);
[m, T] = size(data{1});
K = length(mdl.Theta);
[m_mdl, ncol] = size(mdl.Theta{1});
assert(m == m_mdl, 'Arguments data and mdl are not compatible.')
p = (ncol - 1) / m;

% Model parameters
Tau = mdl.Tau;
Phi = mdl.Phi;

% Compute log-likelihood
logL = 0;
for k = 1:K
    for n = 1:N
        Dnk = dissimMeasure(data, mdl, [n,k]);
        logL = logL + Tau(n,k) * (m/2 * (p-T) * log(2*pi) - ...
                                  .5 * Dnk);
    end
end

% Compute BIC value
score = -2 * logL + (K* ((p + .5) * m^2 + 3/2 * m) + N) * ...
        log(N*(T-p));

end % END of calcBIC


% ================================================================
% Local Functions
% ================================================================
function Dnk = dissimMeasure(data, mdl, nk)
% Compute the dissimilarity measure.

    n = nk(1); k = nk(2);
    N = length(data);
    [m, T] = size(data{1});
    K = length(mdl.Theta);
    [m_mdl, ncol] = size(mdl.Theta{1});
    p = (ncol - 1) / m;

    Dnk = (T-p) * logdet(mdl.Omega{k}) + mdl.Phi(n,k);
end

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
