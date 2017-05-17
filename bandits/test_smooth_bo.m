%% Testing smooth regret for Bayesian Optimization

% Here we compare the regret and switching costs for Thompson
% and cautious thompson (or N-thompson)

clc; clear; close all;

dim = 2;
% average over M experiments
M = 10;
% horizon, i.e. total num of time stages
N = 50;
% GP hyperparameters
hp.type = 'squared exponential ard';
hp.l = 0.2;
hp.scale = 1;
hp.noise.var = 0.01;
ker = @(x1, x2) hp.scale * kernel(x1,x2,hp.l,hp.type);
% buffer for cautious thompson
meshsize = 10;
buffer = floor(0.1*meshsize^2);
delta_gpmi = 10^(-5);

num_algs = 7;
strategy{1}.name = 'GP-UCB';
strategy{1}.delta = 0.1;
strategy{2}.name = 'GP-MI';
strategy{2}.alpha = log(1/delta_gpmi);
strategy{3}.name = 'EI';
strategy{4}.name = 'Thompson-Normal';
strategy{5}.name = 'Thompson-Cautious';
strategy{5}.buffer = buffer;
strategy{6}.name = 'Thompson-Regular';
strategy{6}.lambda = 1/(meshsize^2);
strategy{7}.name = 'Thompson-Plan';
strategy{7}.buffer = buffer;

if dim == 1
    mesh = linspace(0,1.0,meshsize);
elseif dim == 2
    z = linspace(0,1.0,meshsize);
    [X1,X2] = meshgrid(z,z);
    mesh = [X1(:),X2(:)]';
end

regret = zeros(num_algs,N);
cum_regret = zeros(num_algs,N);
idx = zeros(num_algs,1);
switch_cost = zeros(num_algs,N);
cum_switch_cost = zeros(num_algs,N);
x = zeros(dim,num_algs); % actions
last_x = zeros(dim,num_algs); % last actions

for j = 1:M % for each experiment
    
    fprintf('Trial %d\n', j);
    % generating a random function
    
    % initialize different bandit strategies
    for k = 1:num_algs
        bo{k} = BO(mesh,hp,strategy{k});
    end

    % specify variance
    Sigma = ker_matrix(mesh,ker);
    %f = chol(Sigma)*randn(meshsize,1);
    % better for too smooth kernels
    [U,S] = eig(Sigma);
    f = U * sqrt(max(0,S)) * randn(length(mesh),1);
    
    for i = 1:N        
        % play bandit strategy
        for k = 1:num_algs
            [x(:,k),ind] = bo{k}.play();
            y = f(ind) + sqrt(hp.noise.var) * randn;
            % generate noisy function value y
            bo{k}.update(x(:,k),y);
            % calculate regret
            regret(k,i) = max(f) - y;
            switch_cost(k,i) = norm(x(:,k) - last_x(:,k));
            cum_regret(k,i:end) = cum_regret(k,i:end) + regret(k,i);
            cum_switch_cost(k,i:end) = cum_switch_cost(k,i:end) + switch_cost(k,i);
        end
        last_x = x;
    end
end

% plot cumulative regret
figure
subplot(2,1,1);
plot(cum_regret'/M);
legend(strategy{1}.name,...
    strategy{2}.name,...
    strategy{3}.name,...
    strategy{4}.name,...
    strategy{5}.name,...
    strategy{6}.name,...
    'Location','northwest');
title('Growth of cumulative regret');
subplot(2,1,2);
plot((cum_regret + cum_switch_cost)'/M);
legend(strategy{1}.name,...
    strategy{2}.name,...
    strategy{3}.name,...
    strategy{4}.name,...
    strategy{5}.name,...
    strategy{6}.name,...
    'Location','northwest');
title('Growth of cum switching cost + cum regret');
% plot switching costs
