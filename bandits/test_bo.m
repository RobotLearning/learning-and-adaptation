%% Test Bayesian Optimization
% 1 dimensional and no context

clc; clear; close all;

% average over M experiments
M = 20;
% horizon, i.e. total num of time stages
N = 50;
% GP hyperparameters
hp.type = 'squared exponential ard';
hp.l = 0.10;
hp.scale = 1;
hp.noise.var = 0.01;
% hoeffding prob. constant (confidence parameter) for GP-MI
delta_gpmi = 10^(-5);

numAlgs = 4;
strategy{1}.name = 'GP-UCB';
strategy{1}.delta = 0.1;
strategy{2}.name = 'GP-MI';
strategy{2}.alpha = log(1/delta_gpmi);
strategy{3}.name = 'EI';
strategy{4}.name = 'Thompson-Normal';
meshsize = 100;
mesh = linspace(0,1,meshsize);
regret = zeros(numAlgs,N);
cum_regret = zeros(numAlgs,N);
idx = zeros(numAlgs,1);

for j = 1:M % for each experiment
    
    fprintf('Trial %d\n', j);
    % generating a random function
    
    % initialize different bandit strategies
    for k = 1:numAlgs
        bo{k} = bayes_opt(mesh,hp,strategy{k});
    end
    
    % generate function
    % specify variance
    ker = @(x1, x2) hp.scale * kernel(x1,x2,hp.l,hp.type);
    Sigma = ker_matrix(mesh,ker);
    %f = chol(Sigma)*randn(meshsize,1);
    % better for too smooth kernels
    [U,S] = eig(Sigma);
    f = U * sqrt(max(0,S)) * randn(meshsize,1);
    
    for i = 1:N        
        % play bandit strategy
        for k = 1:numAlgs
            [x,ind] = bo{k}.play();
            y = f(ind) + sqrt(hp.noise.var) * randn;
            % generate noisy function value y
            bo{k}.update(x,y);
            % calculate regret
            regret(k,i) = max(f) - y;
            cum_regret(k,i:end) = cum_regret(k,i:end) + regret(k,i);
        end
        
    end
end

% plot cumulative regret
plot(cum_regret'/M);
legend('GP-UCB','GP-MI','EI','Thompson');