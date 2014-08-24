%% Distribution of cumulative regret
% Tries to generate the histogram of cumulative regrets
% 1 dimensional and no context

%% Simulation parameters
clear; clc; close all; 
% number of time steps per run
horizon = 10;
% dimension of the functions to be optimized
dim = 1;
% variance of the noisy function evaluations
var_noise = 0.025;
% generate a fnc over a mesh of pts from that kernel
meshsize = 100;
mesh = linspace(0,1,meshsize);
% number of trials to consider
episode = 5;
% cumulative regret for not picking the maximum
cum_regret = zeros(episode,horizon);
% rewards for each time step, each trial
rewards = zeros(episode,horizon);
% control actions
x = zeros(dim,horizon);

%% Create the function
% lengthscale
l = 0.1;
% scale of the covariance
sigma_s2 = 1;
% sigma of the noisy function evaluations
sigma_noise = sqrt(var_noise);
type = 'squared exponential iso';
kerh = @(x1, x2) sigma_s2 * kernel(x1,x2,l,type);
% specify mean function
mu = zeros(meshsize,1);
% specify variance
Sigma = ker_matrix(mesh,kerh); 
% requires statistical toolbox
fun = mvnrnd(mu, Sigma, episode);
% plot some of the functions
figure(1);
plotsize = 5;
plot(mesh,fun(1:plotsize,:));
names = cell(1,plotsize);
for i = 1:plotsize, names{i} = strcat('function', num2str(i)); end
legend(names);

%% Gaussian Process Optimization

for i = 1:episode
    fprintf('Trial number: %d. \n', i);
    Kinv = 1/(1+var_noise); 
    regret = 0;
    for t = 1:horizon
        % sample x
        [ind,x(t),Kinv] = gpucb(t-1,rewards(i,1:t-1),x(1:t-1),Kinv,kerh);
        % calculate the cost and add noise
        rewards(i,t) = fun(i,ind) + sigma_noise * randn(1);
        % regret for not picking the optimum
        regret = regret + max(fun(i,:)) - fun(i,ind);
        cum_regret(i,t) = regret;                          
    end
end

%% Histogram for cumulative regret

% plot average cumulative regret
figure;
plot(1:horizon,mean(cum_regret));

R_T = cum_regret(:,end);
figure;
nbins = 10;
hist(R_T,nbins);
