% Testing (Multivariate) Linear Bayes Regression

% number of parameters
n = 20;
% number of data points each round (K rounds)
N = 20;

% prior precision matrix
Gamma = 100.0 * eye(n);
beta_nom = randn(n,1);
beta = beta_nom + sqrt(Gamma)\randn(n,1);
var_nom = rand(1);
var = var_nom + rand(1);
mu = beta_nom;
% prior inv gamma distribution values
a = 2; 
b = var_nom;

% draw some values
% make K experiments
K = 10;
err_mean = zeros(K,1);
err_var = zeros(K,1);
log_det_inv_gamma = zeros(K,1);
Xs = [];
ys = [];
for i = 1:K
    X = rand(N,n);
    eps = sqrt(var) * randn(N,1);
    y = X * beta + eps;
    % calculate posterior
    a = a + n/2;
    b = b + 0.5 * (y'*y + mu'*Gamma*mu);
    mu = (X'*X + Gamma)\(Gamma*mu + X'*y);
    Gamma = (X'*X + Gamma);
    b = b - 0.5 * (mu'*Gamma*mu);
    % calculate error
    err_mean(i) = norm(beta - mu,2);
    err_var(i) = norm((b/(a-1)) - var,2);
    % how big is the new mean estimate variance
    log_det_inv_gamma(i) = log(1/det(Gamma));
    % save for batch processing later
    Xs = [Xs; X];
    ys = [ys; y];
end

% batch regression for comparison
mu_batch = Xs \ ys;
err_batch = norm(mu_batch - beta,2);