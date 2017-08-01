%% Testing smooth regret for Bayesian Optimization

% Here we compare the regret and switching costs for Thompson
% and cautious thompson (or N-thompson)

clc; clear; close all;

rng(3);
dim = 2;
% average over M experiments
M = 5;
% horizon, i.e. total num of time stages
T = 30;
% discretization
N = 10;

if dim == 1
    mesh = linspace(0,1.0,N);
elseif dim == 2
    z = linspace(0,1.0,N);
    [X1,X2] = meshgrid(z,z);
    mesh = [X1(:),X2(:)]';
end

% GP hyperparameters
hp.type = 'squared exponential ard';
hp.l = 0.2;
hp.scale = 1;
hp.noise.var = 0.1;
ker = @(x1, x2) hp.scale * kernel(x1,x2,hp.l,hp.type);
Sigma = ker_matrix(mesh,ker);
[U,S] = eig(Sigma);
% buffer for cautious thompson
lambda = 5;
buffer = @(t) ceil(length(mesh)/(10*t));
delta_gpmi = 10^(-5);

% Mismatch
% hp2 = hp;
% hp2.l = 0.15;

num_algs = 6;
strategy{1}.name = 'Thompson-Cautious';
strategy{1}.buffer = buffer;
strategy{2}.name = 'Thompson-Normal';
strategy{3}.name = 'Thompson-Regular';
strategy{3}.lambda = lambda/length(mesh);
strategy{4}.name = 'GP-UCB';
strategy{4}.delta = 0.1;
strategy{5}.name = 'GP-MI';
strategy{5}.alpha = log(1/delta_gpmi);
strategy{6}.name = 'EI';
% strategy{7}.name = 'Thompson-Plan';
% strategy{7}.buffer = buffer;

regret = zeros(num_algs,T);
cum_regret = zeros(num_algs,T);
idx = zeros(num_algs,1);
switch_cost = zeros(num_algs,T);
cum_switch_cost = zeros(num_algs,T);
x = zeros(dim,num_algs); % actions
last_x = zeros(dim,num_algs); % last actions

for j = 1:M % for each experiment
    
    fprintf('Trial %d\n', j);
    % generating a random function
    
    % initialize different bandit strategies
    for k = 1:num_algs
        bo{k} = BO(mesh,hp,strategy{k});
    end

    % SAMPLING FROM GP
    %f = U * sqrt(max(0,S)) * randn(length(mesh),1);
    % USING BRANIN FUNCTION
    f = branin(mesh);
    f = f/min(f);
    
    for i = 1:T
        % play bandit strategy
        for k = 1:num_algs
            [x(:,k),ind] = bo{k}.play();
            y = f(ind) + sqrt(hp.noise.var) * randn;
            % generate noisy function value y
            bo{k}.update(x(:,k),y);
            % calculate regret
            regret(k,i) = max(f) - y;
            switch_cost(k,i) = lambda * norm(x(:,k) - last_x(:,k));
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
title('Average cumulative regret');
subplot(2,1,2);
plot((cum_regret + cum_switch_cost)'/M);
legend(strategy{1}.name,...
    strategy{2}.name,...
    strategy{3}.name,...
    strategy{4}.name,...
    strategy{5}.name,...
    strategy{6}.name,...
    'Location','northwest');
title('Average cumulative augmented regret');
% plot switching costs
