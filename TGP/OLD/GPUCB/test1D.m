%% 1-D test case
% Added mismatch comparison
% Adapted from Andreas Krause's paper

clear; clc; close all; 

%%%%%%%%%%%%%%%%%%%% SIMULATION PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of time steps per run
horizon = 50;
% dimension of the functions to be optimized
dim = 1;
% variance of the noisy function evaluations
var_noise = 0.025;
% generate a fnc over a mesh of pts from that kernel
meshsize = 100;
mesh = linspace(0,1,meshsize);
% averaging over number of trials
trials = 5;
% cumulative regret for not picking the optimum
cum_regret = zeros(trials,horizon);
% costs for picking the control action
cost = zeros(trials, horizon);
% control actions
x = zeros(dim,horizon);

%% FUNCTION CREATION

%%%%%%%%%%%%%%%%%%%%%%%%%% SPECIFY KERNEL STRUCTURE %%%%%%%%%%%%%%%%%%%%%%%
% lengthscale
l = 0.1;
% scale of the covariance
sigma_s2 = 1;
% sigma of the noisy function evaluations
sigma_noise = sqrt(var_noise);
type = 'squared exponential iso';
ker_handle = @(x1, x2) sigma_s2 * kernel(x1,x2,l,type);

%%%%%%%%%%%%%%%%%%%%%%%%%% CONSIDER MISMATCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% control actions under mismatch
x_mm = zeros(dim,horizon);
% costs under mismatch
cost_mm = zeros(trials, horizon);
% cumulative regret for not picking the optimum under mismatch
cum_regret_mm = zeros(trials,horizon);
% lengthscale
l_mm = 0.1;
sigma_s2_mm = 2;
eps_mm = 0; %0.025;
var_noise_mm = var_noise + eps_mm;
ker_handle_mm = @(x1,x2) sigma_s2_mm * kernel(x1,x2,l_mm,type);


%%%%%%%%%%%%%%%%%%%%%%%%%% GENERATE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specify mean function
mu = zeros(meshsize,1);
% specify variance
Sigma = ker_matrix(mesh,ker_handle); 
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

% plot some of the functions
figure(1);
plotsize = 5;
plot(mesh,fun(1:plotsize,:));
names = cell(1,plotsize);
for i = 1:plotsize, names{i} = strcat('function', num2str(i)); end
legend(names);

%%%%%%%%%%%%%%%%%%%% PASS TO STRUCTURE %%%%%%%%%%%%%%%%%
% create structure to pass parameters easily
STR.dim = 1;
STR.meshsize = meshsize;
STR.lb = min(mesh);
STR.ub = max(mesh);
STR.l = l; % kernel parameter needed in derivative

%% GP-UCB 
tic;
for i = 1:trials
    fprintf('Trial number: %d. \n', i);
    K = 0; 
    K_mm = 0; % parameter mismatch
    %Kinv = 1/(1+var_noise); 
    %Kinv_mm = 1/(1+var_noise_mm);
    regret = 0;
    regret_mm = 0;
    for t = 1:horizon
        % compute control input u
        STR.kernel = ker_handle;
        STR.var = var_noise;
        STR.K = K;
        [x(t), STR] = gp_ucb(t-1, cost(i,1:t-1), x(1:t-1), STR);
        K = STR.K;
        % calculate the cost and add noise
        cost(i,t) = interp1(mesh, fun(i,:), x(t)) + sigma_noise * randn(1);
        % regret for not picking the optimum
        regret = regret + interp1(mesh, fun(i,:), x(t)) - min(fun(i,:));
        cum_regret(i,t) = regret;                          
        
        % PARAMETER MISMATCH
        %%{
        % compute control input u under hyperparameter mismatch
        STR.var = var_noise_mm;
        STR.kernel = ker_handle_mm;
        STR.K = K_mm;
        [x_mm(t), STR] = gp_ucb(t-1, cost_mm(i,1:t-1), x_mm(1:t-1), STR);
        K_mm = STR.K;
        % calculate the cost and add noise
        cost_mm(i,t) = interp1(mesh, fun(i,:), x_mm(t)) + sigma_noise * randn(1);
        % regret for not picking the optimum
        regret_mm = regret_mm + interp1(mesh, fun(i,:), x_mm(t)) - min(fun(i,:));
        cum_regret_mm(i,t) = regret_mm;              
        % debug : inform the user when it's slow
        %str = fprintf('At time stage %d, control: %f, cumulative regret: %f.\n', ...
        %    t, x(t), cum_regret(i,t));
        %}
    end
end

%%%%%%%%%%%%%%%%%%%%%%% INFORM USER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('GP-UCB took %f seconds to run on %d iterations with %d horizon.\n', ...
        toc, trials, horizon);
%plot cumulative regret
Ravg = sum(cum_regret,1)/trials;
Ravg_mm = sum(cum_regret_mm,1)/trials;
figure(2);
plot(1:horizon, Ravg, '-r', 1:horizon, Ravg_mm);
legend('Regret under zero mismatch', 'Regret under mismatch');