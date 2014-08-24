function [hyp,bhat] = reml(hyp, covfunc, likfunc, x, y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A very crude implementation of REML 
% Uses gp and minimize functions in Rasmussen's toolbox 
% Feeds an additional complexity/derivative complexity term to the 
% conjugate gradient (CG) minimization at each step n of the EM algorithm
%
% Author : Okan Koc, Msc ETH Zurich
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% transform data
Xbar = [ones(length(x),1), x]; % design matrix

num_theta = length(hyp.cov) + 1; % sigma_n also included
num_beta = size(Xbar,2);
change = Inf(num_theta + num_beta,1);
eps = 1e-2;
iter = 0;
bhat_old = hyp.mean;
hyp.mean = [];
Ncg = 1000;
while max(change) > eps
    % covariance matrix 
    K = feval(covfunc{:}, hyp.cov, x) + ...
        (exp(hyp.lik)^2)*eye(length(x));
    pseudo = Xbar' * (K \ Xbar);
    bhat = pseudo \ (Xbar' * (K \ y)); %mean hps
    % transformed y-values
    y_tr = y - Xbar*bhat;
    % save old value
    hyp_old = hyp;
    % do Full Maximum Likelihood Estimation on y_tr
    hyp = minimin(hyp, @gpMod, -Ncg, false, @infExact, [], covfunc, ...
               likfunc, x, y_tr);
    change(1:length(hyp.cov)) = abs(exp(hyp.cov) - exp(hyp_old.cov));
    change(length(hyp.cov)+1) = abs(exp(hyp.lik) - exp(hyp_old.lik));
    change(length(hyp.cov)+2:end) = abs(bhat - bhat_old);
    bhat_old = bhat;
    iter = iter + 1;
end
fprintf('REML converged in %d iterations. \n', iter);

end

function [nll, dnll] = gpMod(hyp, inf, mean, cov, lik, x, y)

[nll_ML, dnll_ML] = gp(hyp, inf, mean, cov, lik, x, y);

Xbar = [ones(length(x),1), x]; % design matrix
K = feval(cov{:}, hyp.cov, x);
Knoisy = K + (exp(hyp.lik)^2)*eye(length(x));
% complexity penalty
mat = Knoisy \ Xbar;
pseudo = Xbar' * mat;
comp = log(det(pseudo))/2;
% compute complexity derivative
dKn = zeros(length(K), length(K), length(hyp.cov)+1,1); 
diffX = zeros(length(K),length(K),size(x,2));
% derivative for length scales
for d = 1:size(x,2)
    % form multiplicative matrix
    for j = 1:size(x,1)
        diffX(:,j,d) = (x(:,d) - x(j,d)).^2;
    end
    dKn(:,:,d) = (1/(exp(hyp.cov(d))^3))*diffX(:,:,d).*K;
end
% for sigma_s
dKn(:,:,length(hyp.cov)) = (2/exp(hyp.cov(end)))*K;
% for sigma n
dKn(:,:,end) = 2*exp(hyp.lik)*eye(length(K));

% use dKn to compute complexity derivative
comp_dnll.mean = [];
comp_dnll.cov = zeros(length(hyp.cov),1);
comp_dnll.lik = 0;
for i = 1:length(hyp.cov)
    comp_dnll.cov(i) = -sum(diag((pseudo\(mat'*dKn(:,:,i)*mat))));
end
comp_dnll.lik = -sum(diag(pseudo\(mat'*dKn(:,:,end)*mat)));

% output: value and derivative of negative log likelihood
nll = nll_ML + comp;
dnll.mean = dnll_ML.mean + comp_dnll.mean/2;
dnll.cov = dnll_ML.cov + comp_dnll.cov/2;
dnll.lik = dnll_ML.lik + comp_dnll.lik/2; 

end