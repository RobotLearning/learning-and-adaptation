%% Mean and Distribution of cumulative regret
% for Gaussian Process Optimization
% 1 dimensional and no context

%% Simulation parameters
clear; clc; close all; 
% number of time steps per run
horizon = 50;
% dimension of the functions to be optimized
dim = 1;
% variance of the noisy function evaluations
var_noise = 0.025;
% generate a fnc over a mesh of pts from that kernel
meshsize = 100;
mesh = linspace(0,1,meshsize);
% number of trials to consider
trials = 5;
% rewards for each time step, each trial
rewards = zeros(trials,horizon);
% control actions
x = zeros(dim,horizon);
% get the same results always
rng(0);

%% Create the functions
% lengthscale
l = 0.1;
% scale of the covariance
sigma_s2 = 1;
% sigma of the noisy function evaluations
sigma_noise = sqrt(var_noise);
type = 'squared exponential iso';
ker = @(x1, x2) sigma_s2 * kernel(x1,x2,l,type);
% specify mean function
mu = zeros(meshsize,1);
% specify variance
Sigma = ker_matrix(mesh,ker); 
% requires statistical toolbox
if exist('mvnrnd')
    fun = mvnrnd(mu, Sigma, trials);
else
    for i = 1:trials
        lambda = 1e-3;
        fun(:,i) = mu + chol(Sigma + lambda*eye(meshsize))*randn(meshsize,1);
    end
    fun = fun';
end
% plot some of the functions
figure(1);
plotsize = min(5,trials);
plot(mesh,fun(1:plotsize,:));
names = cell(1,plotsize);
for i = 1:plotsize, names{i} = strcat('function', num2str(i)); end
legend(names);

%% Bandit Optimization using GPUCB criterion

% cumulative regret for not picking the maximum
ucb_regret = zeros(trials,horizon);
% CONSTRUCT BETA CONSTANT 
% hoeffding prob. constant
delta = 0.1;

% GPUCB Thm1 for discrete space
beta_thm1 = @(t,D,delta) 2*log((t^2)*D*(pi^2)/(6*delta));
disp('Starting GPUCB experiment...');

for i = 1:trials
    fprintf('Trial number: %d. \n', i);
    Kinv = 1/(1+var_noise); 
    regret = 0;
    for t = 1:horizon
        % construct acquisition function
        beta = beta_thm1(t,meshsize,delta);
        acq = @(m,s2) m + sqrt(beta * s2);
        % sample x
        [ind,x(t),Kinv,~] = bandit(t-1,rewards(i,1:t-1),x(1:t-1),Kinv,...
                                 ker,mesh,acq);
        % calculate the cost and add noise
        rewards(i,t) = fun(i,ind) + sigma_noise * randn(1);
        % regret for not picking the optimum
        regret = regret + max(fun(i,:)) - fun(i,ind);
        ucb_regret(i,t) = regret;                          
    end
end

%% Bandit Optimization using MI criterion

% cumulative regret for not picking the maximum
mi_regret = zeros(trials,horizon);
% rewards for each time step, each trial
rewards = zeros(trials,horizon);
% control actions
x = zeros(dim,horizon);

% hoeffding prob. constant (confidence parameter)
delta = 10^(-9);
% alpha parameter for MI acquisition function
alpha = log(1/delta);

disp('Starting GP-MI experiment...');

for i = 1:trials
    fprintf('Trial number: %d. \n', i);
    Kinv = 1/(1+var_noise); 
    regret = 0;
    gamma = 0;
    for t = 1:horizon
        % construct acquisition function
        acq = @(m,s2) m + sqrt(alpha * (s2 + gamma)) - sqrt(alpha * gamma);
        % sample x
        [ind,x(t),Kinv,s2max] = bandit(t-1,rewards(i,1:t-1),x(1:t-1),Kinv,...
                                 ker,mesh,acq);
        % calculate the cost and add noise
        rewards(i,t) = fun(i,ind) + sigma_noise * randn(1);
        % regret for not picking the optimum
        regret = regret + max(fun(i,:)) - fun(i,ind);
        mi_regret(i,t) = regret;
        % update gamma
        gamma = gamma + s2max;
    end
end

%% Bandit Optimization using EI criterion

% cumulative regret for not picking the maximum
ei_regret = zeros(trials,horizon);
% rewards for each time step, each trial
rewards = zeros(trials,horizon);
% control actions
x = zeros(dim,horizon);

disp('Starting Expected Improvement experiment...');
density = @(x) 1/sqrt(2*pi) .* exp((-1/2).*(x.^2));
distr = @(x) (1 - erf(x./sqrt(2)))./2;

for i = 1:trials
    fprintf('Trial number: %d. \n', i);
    Kinv = 1/(1+var_noise); 
    regret = 0;
    for t = 1:horizon
        % construct acquisition function
        acq = @(m,s2) (m - max(m)) .* distr((max(m) - m)./sqrt(s2)) + ...
                      sqrt(s2) .* density((max(m) - m)./sqrt(s2));
        % sample x
        [ind,x(t),Kinv,~] = bandit(t-1,rewards(i,1:t-1),x(1:t-1),Kinv,...
                                 ker,mesh,acq);
        % calculate the cost and add noise
        rewards(i,t) = fun(i,ind) + sigma_noise * randn(1);
        % regret for not picking the optimum
        regret = regret + max(fun(i,:)) - fun(i,ind);
        ei_regret(i,t) = regret;                  
    end
end

%% Histogram for cumulative regret

% plot average cumulative regret
figure(2);
plot(1:horizon,mean(ucb_regret),...
     1:horizon,mean(ei_regret),...
     1:horizon,mean(mi_regret));
title('Average cumulative regret over time');
legend('GPUCB','EI','GP-MI');

% R_T = cum_regret(:,end);
% figure(3);
% nbins = 10;
% hist(R_T,nbins);
% title('Histogram of cumulative regrets in 10 bins');
