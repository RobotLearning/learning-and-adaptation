%% Multi-D test case
% Adapted from Andreas Krause's paper

clear; clc; close all; 

%%%%%%%%%%%%%%%%%%%% SIMULATION PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of time steps per run
horizon = 100;
% dimension of the functions to be optimized
dim = 3;
% variance of the noisy function evaluations
var_noise = 0.01;
% averaging over number of trials
trials = 5;
% cumulative regret for not picking the optimum
cum_regret = zeros(trials,horizon);
% costs for picking the control action
cost = zeros(trials, horizon);
% control actions
x = zeros(dim,horizon);

%%%%%%%%%%%%%%%%%%%%%%%%% GENERATE MESH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate a fnc over a mesh of pts from that kernel
n = 7; % in one dimension
meshsize = n^dim; 
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
Lambda = l1 * (1:dim); % increasing regularity in dimensions
% scale of the covariance
sigma_s2 = 1;
% sigma of the noisy function evaluations
sigma_noise = sqrt(var_noise);
type = 'squared exponential ard';
ker_handle = @(x1, x2) sigma_s2 * kernel(x1,x2,Lambda,type);

%%%%%%%%%%%%%%%%%%%%%%%%%% GENERATE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specify mean function
mu = zeros(meshsize,1);
% specify variance
Sigma = ker_matrix(X',ker_handle); 
% evaluate the function

% FIRST METHOD - Sigma must be positive definite
%{
fun = zeros(trials, meshsize);
% small epsilon to keep matrix positive definite
eps = 0.003;
for i = 1:trials
    fun(i,:) = mu(:) + chol(Sigma + eps*eye(meshsize)) * randn(meshsize,1);
end
%}
% SECOND METHOD - requires statistical toolbox though!
fun = mvnrnd(mu, Sigma, trials);

% plot the first function generated
if dim == 2
    figure(1);
    Z = zeros(n,n,trials);
    for i = 1:trials
        Z(:,:,i) = reshape(fun(i,:), n, n);
    end
    surf(XX{1},XX{2},Z(:,:,1));
end

% arrange Z array
ns = cell(1,dim);
Z = cell(1,trials);
[ns{:}] = deal(n);
for i = 1:trials
    Z{i} = reshape(fun(i,:), ns{:});
end    

%%%%%%%%%%%%%%%%%%%% PASS TO STRUCTURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create structure to pass parameters easily
STR.dim = dim;
STR.meshsize = meshsize;
STR.lb = min(xx) * ones(dim,1);
STR.ub = max(xx) * ones(dim,1);
STR.l = Lambda; % kernel parameter needed in derivative
STR.sigma_s2 = sigma_s2;

%% GP-UCB 
inp = cell(1,dim);
tic;
for i = 1:trials
    fprintf('Trial number: %d. \n', i);
    K = 0; 
    %Kinv = 1/(1+var_noise); 
    regret = 0;
    for t = 1:horizon
        % compute control input u
        STR.kernel = ker_handle;
        STR.var = var_noise;
        STR.K = K;
        [x(:,t), STR] = gp_ucb(t-1, cost(i,1:t-1), x(:,1:t-1), STR);
        K = STR.K;
        for j = 1:dim, inp{j} = x(j,t); end
        % calculate the cost and add noise
        cost(i,t) = interpn(XX{:}, Z{i}, inp{:}) + sigma_noise * randn(1);
        % regret for not picking the optimum
        regret = regret + interpn(XX{:}, Z{i}, inp{:}) - min(fun(i,:));
        cum_regret(i,t) = regret;                                 
    end
end

%%%%%%%%%%%%%%%%%%%%%%% INFORM USER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('GP-UCB took %f seconds to run on %d iterations with %d horizon.\n', ...
        toc, trials, horizon);
%plot cumulative regret
Ravg = sum(cum_regret,1)/trials;
figure(2);
plot(1:horizon, Ravg);
legend('Regret under zero mismatch');