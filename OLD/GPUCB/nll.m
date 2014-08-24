function [val, der] = nll(theta, x, y, type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Negative Log Likelihood given parameter values
%
% Inputs : x - data points, N x d matrix, d : dimensions, N : sample size
%          y - observations, N x m matrix with m different experiments
%          theta - hyperparameter values
%          type - type of kernel
%
% Outputs: val - negative log likelihood
%          der - derivative of nll
%
% Author : Okan Koc, IDSC ETH Zurich
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% number of dimensions
d = size(x,2);
% number of data points
N = size(x,1);
% number of experiments
m = size(y,2);

% extract hyperparameters
Lambda = theta(1:d);
sigma_s = theta(d+1);
sigma_s2 = sigma_s^2;
sigma_n = theta(end);
sigma_n2 = sigma_n^2;
ker_handle = @(x1, x2) sigma_s2 * kernel(x1,x2,Lambda,type);

% construct covariance matrix
K = ker_matrix(x',ker_handle); 
Kn = K + sigma_n2 * eye(length(K));

% negative log likelihood
% assuming mu is zero
mat = Kn \ y;
% numerically stable way to compute log of det
logdetKn = sum(log(eig(Kn)));
val = m*logdetKn/2 + sum(diag(y' * mat)); 

% derivative of negative log likelihood

% compute derivative of Kn
dKn = zeros(length(K), length(K), d+2); 
diffX = zeros(length(K),length(K),d);
% derivative for length scales
for i = 1:d
    % form multiplicative matrix
    for j = 1:N
        diffX(:,j,d) = (x(:,d) - x(j,d)).^2;
    end
    dKn(:,:,d) = (1/Lambda(i)^3)*diffX(:,:,d).*K;
end
% for sigma_s
dKn(:,:,d+1) = (2/sigma_s)*K;
% for sigma n
dKn(:,:,end) = 2*sigma_n*eye(length(K));
    
der = zeros(d+2,1);
for i = 1:d+2
    der(i) = m*sum(diag(Kn \ dKn(:,:,i)))/2 - ...
             sum(diag(mat' * dKn(:,:,i) * mat));
end