%% Checking consistency of Maximum Likelihood Estimation

clear; clc; close all; 

%%%%%%%%%%%%%%%%%%%% SIMULATION PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dimension of the functions to be optimized
dim = 1;
% variance of the noisy function evaluations
var_noise = 0.01;
% number of training sets for likelihood estimation
train = 1;

%%%%%%%%%%%%%%%%%%%%%%%%% GENERATE MESH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate a fnc over a mesh of pts from that kernel
n = 50; 
meshsize = n^dim; % in one dimension
xx = linspace(0,1,n);
% generate n ^dim training points
[varargout{1:dim}] = ndgrid(xx);
XX = varargout;
X = zeros(meshsize,dim);
% vectorize cells
for i = 1:dim
    X(:,i) = XX{i}(:);
end

%% FUNCTION CREATION

%%%%%%%%%%%%%%%%%%%%%%%%%% SPECIFY KERNEL STRUCTURE %%%%%%%%%%%%%%%%%%%%%%%
% lengthscales
l1 = 0.1;
l2 = 0.2;
%Lambda = [l1,l2];
Lambda = l1 * (1:dim);
% scale of the covariance
sigma_s2 = 1;
% sigma of the noisy function evaluations
sigma_noise = sqrt(var_noise);
type = 'squared exponential ard';
ker_handle = @(x1, x2) sigma_s2 * kernel(x1,x2,Lambda,type);


%% MAXIMUM LIKELIHOOD ESTIMATION (ML)

% specify mean function
mu = zeros(1,meshsize);
% specify variance
Sigma = ker_matrix(X',ker_handle); 
% evaluate the function
%fun_train = zeros(meshsize, train);
% small epsilon to keep matrix positive definite
eps = var_noise;
%for i = 1:train
%    fun_train(:,i) = mu(:) + chol(Sigma + eps*eye(meshsize)) * randn(meshsize,1);
%end
fun_train = mvnrnd(mu, Sigma, train) + sigma_noise*randn(train,meshsize);

% Derivative Based Optimization
% run conjugate gradient method to find best set of parameters
disp('MLE: Derivative based optimization...');
p0 = 0.1;
hpm = p0 * ones(dim+2,1);
Ncg = 100;
hpm = minimin(hpm, @nll, -Ncg, true, X, fun_train', type);
hpm = abs(hpm);

% print inferred values
str = [];
ls = cell(1,dim);
for i = 1:dim, 
    str = [str, ['ell', num2str(i), ' = %f, ']]; %#ok
    ls{i} = hpm(i);
end 
str = [str, 'sigma_s = %f, sigma_n = %f \n'];
fprintf(str, ls{:}, hpm(end-1), hpm(end));

% plot nll of different l1 values
%{
func = @(hpm) nll(hpm,X',fun_train',type);
l1vals = linspace(0.01,1,50);
[nlls,dnlls] = arrayfun(func,l1vals);
% numerical derivative
dnll_num = diff(nlls)./diff(l1vals);
dnll_num(length(l1vals)) = dnll_num(end);
plot(l1vals,nlls,'-r',l1vals,dnlls,'-b',l1vals,dnll_num,'-g');
legend('nll values', 'nll derivatives', 'numerical derivatives');
%}

%{
% Sampling Based Optimization
disp('MLE: Sampling based optimization...');
% Using DIRECT 
% Send options to Direct 
options.showits   = 0;
options.tol       = 0.05;
options.maxevals  = 100;
options.maxits    = 50;
options.maxdeep   = 100;
% Pass function as part of a Matlab Structure

lb = l1vals(1);
ub = l1vals(end);
bounds = [lb, ub];
problem.f = func;
[~,hpm,hist] = Direct(problem,bounds,options);
fprintf('ell = %f \n', hpm(1));
%}