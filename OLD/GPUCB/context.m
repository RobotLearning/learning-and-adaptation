%% - CONTEXTUAL BANDITS
% Adapted from Andreas Krause's paper
%
% TODOS:
%
% 5 - Try known smooth functions (brainin, sines, ...)
% 5.2 - RKHS approach ??

clear; clc; close all; 

%% SIMULATION PARAMETERS 
% number of time steps per run
horizon = 100;
% dimension of the functions to be optimized
dim_ctx = 1;
dim_act = 2;
dim = dim_ctx + dim_act;
% variance of the noisy function evaluations
var_noise = 0.01;
% averaging over number of trials
trials = 1;
% cumulative regret for not picking the optimum
cum_regret = zeros(trials,horizon);
% costs for picking the control action
cost = zeros(trials, horizon);
% control actions
u = zeros(dim_act,horizon);
% contexts + actions = cxu
cxu = zeros(dim, horizon);

%% GENERATE MESH 
% generate a fnc over a mesh of pts from that kernel
n_ctx = 10; % number of different contexts
n_act = 5; % different actions
meshsize_ctx = n_ctx^(dim_ctx);
meshsize_act = n_act^(dim_act);
meshsize = meshsize_ctx * meshsize_act; 
axis_ctx = linspace(0,1,n_ctx);
axis_act = linspace(0,1,n_act);
xxs_ctx = cell(1,dim_ctx);
xxs_act = cell(1,dim_act);
[xxs_ctx{:}] = deal(axis_ctx);
[xxs_act{:}] = deal(axis_act);
% generate n ^dim training points
[varargout{1:dim}] = ndgrid(xxs_ctx{:}, xxs_act{:});
XX = varargout;
X = zeros(meshsize,dim);
% vectorize cells
for i = 1:dim
    X(:,i) = XX{i}(:);
end

%% SPECIFY KERNELS

%%%%%%%%%%%%%%%%%%%%%%%%%% SPECIFY KERNEL STRUCTURE %%%%%%%%%%%%%%%%%%%%%%%
% lengthscales
l_ctx = 0.2 * (1:dim_ctx); % increasing regularity in dimensions
l_act = 0.1 * (1:dim_act);
Lambda = [l_ctx, l_act];
% scale of the covariance
sigma_s2 = 1;
% sigma of the noisy function evaluations
sigma_noise = sqrt(var_noise);
type_ctx = 'squared exponential ard';
ker_ctx = @(x1,x2) sigma_s2 * kernel(x1,x2,l_ctx,type_ctx);
type_act = 'linear ard';
ker_act = @(x1,x2) kernel(x1,x2,l_act,type_act);
ker_handle = @(x1,x2) ker_ctx(x1(1:dim_ctx),x2(1:dim_ctx))*...
                      ker_act(x1(dim_ctx+1:end),x2(dim_ctx+1:end));
                  
%%%%%%%%%%%%%%%%%%%%%% CONSIDER NOCONTEXT CASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here mismatch is no use of context and using only one GP-UCB
% control actions under mismatch
u_noctx = zeros(dim_act,horizon);
% costs under mismatch
cost_noctx = zeros(trials, horizon);
% cumulative regret for not picking the optimum under mismatch
cum_regret_noctx = zeros(trials,horizon);
ker_handle_noctx = @(x1,x2) sigma_s2 * kernel(x1,x2,l_act,type_act);


%% GENERATE FUNCTIONS

% specify quadratic function for the mean
%%{
m = 10;
Q_boost = m * eye(dim_act);
c_boost = -m * ones(dim_act,1);
a = m/2;
boost = @(X) diag(X'*Q_boost*X) + X'*c_boost + a;
mu = boost(X(:,dim_ctx+1:end)'); mu = mu(:);
%}
% specify zero mean
%{ 
boost = @(X) 0;
Q_boost = 0;
c_boost = 0;
% specify mean function
mu = zeros(meshsize,1);
%}
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
rng(1);
fun = mvnrnd(mu, Sigma, trials);

%% VISUALIZE FIRST FUNCTION
if dim == 2
    figure(1);
    Z = zeros(n_ctx,n_act,trials);
    for i = 1:trials
        Z(:,:,i) = reshape(fun(i,:), n_ctx, n_act);
    end
    surf(XX{1},XX{2},Z(:,:,1));
    xlabel('contexts');
    ylabel('actions');
    zlabel('costs');
end
% visualize mu
% MU = reshape(mu, n_ctx, n_act);
% surf(XX{1}, XX{2}, MU);

% arrange Z array
ns_ctx = cell(1,dim_ctx);
ns_act = cell(1,dim_act);
Z = cell(1,trials);
[ns_ctx{:}] = deal(n_ctx);
[ns_act{:}] = deal(n_act);
for i = 1:trials
    Z{i} = reshape(fun(i,:), ns_ctx{:}, ns_act{:});
end    

%% PASS TO STRUCTURE 
% create structure to pass parameters easily
STR.ker_ctx = ker_ctx;
STR.kernel = ker_handle;
STR.var = var_noise;
STR.meshsize = meshsize;
STR.dim = dim;
STR.dim_ctx = dim_ctx;
STR.dim_act = dim_act;
STR.boost.fun = boost;
STR.boost.Q = Q_boost;
STR.boost.c = c_boost;
STR.l = Lambda; % kernel parameter needed in derivative
STR.l_act = l_act; 
STR.l_ctx = l_ctx;
STR.lb = min(axis_act) * ones(dim_act,1);
STR.ub = max(axis_act) * ones(dim_act,1);
STR.sigma_s2 = sigma_s2;

%% GP-UCB 
inp = cell(1,dim);
tic;
for i = 1:trials
    fprintf('Trial number: %d. \n', i);
    K = 0; 
    ctx = zeros(dim_ctx,horizon);
    %K_noctx = 0;
    Kinv = 1/(1+var_noise); 
    regret = 0;
    %regret_noctx = 0;
    for t = 1:horizon
        % reveal context out of n_ctx points
        ctx_ind = randi(n_ctx,1,dim_ctx);
        ctx(:,t) = axis_ctx(ctx_ind);
        STR.ctx = ctx(:,t); 
        STR.Kinv = Kinv;
        STR.K = K;
        %STR.kernel = ker_handle;
        %STR.var = var_noise;
        %STR.meshsize = meshsize;
        %STR.dim = dim;
        %STR.l = Lambda; % kernel parameter needed in derivative
        [x, STR] = cgp_ucb(t-1, cost(i,1:t-1), ctx(:,1:t-1), ...
                                u(:,1:t-1), STR);
        Kinv = STR.Kinv;
        K = STR.K;
        u(:,t) = x;
        for j = 1:dim_ctx, inp{j} = ctx(j,t); end
        for j = 1:dim_act, inp{dim_ctx+j} = u(j,t); end
        % update past contexts + actions
        %cxu(:,t) = [ctx(:,t); u(:,t)];
        % calculate the cost and add noise
        cost(i,t) = interpn(XX{:}, Z{i}, inp{:}) - boost(u(:,t)) + ...
                    sigma_noise * randn(1);
        % find the minimum at the context
        minval = min(Z{i}(ctx_ind,:));
        % regret for not picking the optimum
        regret = regret + interpn(XX{:}, Z{i}, inp{:}) - minval;
        cum_regret(i,t) = regret;
        
        % display every 10th iteration if process is slow
        if(rem(t,10) == 0), fprintf('Time stage %d ...\n', t); end
       
        % NO CONTEXT 
        %{
        % compute control input u under hyperparameter mismatch
        STR.var = var_noise;
        STR.kernel = ker_handle_noctx;
        STR.K = K_noctx;
        STR.meshsize = meshsize_act;
        STR.dim = dim_act;
        STR.l = l_act;
        [u_noctx(:,t), STR] = gp_ucb(t-1, cost_noctx(i,1:t-1), ...
                                   u_noctx(1:t-1), STR);
        K_noctx = STR.K;
        for j = 1:dim_ctx, inp{j} = ctx(j); end
        for j = 1:dim_act, inp{dim_ctx+j} = u_noctx(j,t); end
        % calculate the cost and add noise
        cost_noctx(i,t) = interpn(XX{:}, Z{i}, inp{:}) + ...
                          sigma_noise * randn(1);
        % regret for not picking the optimum
        regret_noctx = regret_noctx + interpn(XX{:}, Z{i}, inp{:}) ...
                       - minval;
        cum_regret_noctx(i,t) = regret_noctx;              
        %}
    end
end

%% INFORM USER 
fprintf('GP-UCB took %f seconds to run on %d iterations with %d horizon.\n', ...
        toc, trials, horizon);
%plot cumulative regret
Ravg = sum(cum_regret,1)/trials;
%Ravg_noctx = sum(cum_regret_noctx,1)/trials;
figure(2);
plot(1:horizon, Ravg);
legend('CGP-UCB Regret');
%plot(1:horizon, Ravg, '-r', 1:horizon, Ravg_noctx);
%legend('CGP-UCB Regret under zero mismatch', ...
%       'GP-UCB Regret ignoring context');