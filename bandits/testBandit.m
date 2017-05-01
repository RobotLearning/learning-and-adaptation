%% Testing different bandit algorithms 
% 1d and no context

clc; clear; close all;

% average over M trials
M = 100;
% number of arms
K = 50;
% horizon, i.e. total num of time stages
N = 1000;

numAlgs = 4;
strategy{1}.name = 'EPS-GREEDY';
strategy{1}.eps = @(t) 10/t;
strategy{2}.name = 'UCB1';
strategy{2}.rho = 1;
strategy{3}.name = 'UCB1-V';
strategy{4}.name = 'ThompsonGauss';
strategy{4}.a = 10;
strategy{4}.b = 0.5;

regret = zeros(numAlgs,N);
cum_regret = zeros(numAlgs,N);
idx = zeros(numAlgs,1);

for j = 1:M
    
    % generating K gaussians
    mu = 1*rand(K,1);
    var = 1*rand(K,1);
    
    % initialize different bandit strategies
    for k = 1:numAlgs
        bandit{k} = Bandit(K,strategy{k});
    end
    
    for i = 1:N
        % generate rewards
        rewards = mu + sqrt(var).*randn(K,1);
        % generate from uniform distribution
        %rewards = (mu - sqrt(12*var)/2) + sqrt(12*var).*rand(K,1);
        
        % play bandit strategy
        for k = 1:numAlgs
            idx(k) = bandit{k}.play();
            bandit{k}.reward(rewards(idx(k)),idx(k));
            % calculate regret
            regret(k,i) = max(mu) - rewards(idx(k));
            cum_regret(k,i:end) = cum_regret(k,i:end) + regret(k,i);
        end
        
    end
end

% plot cumulative regret
plot(cum_regret'/M);
legend('\epsilon-Greedy','UCB_1','UCB_1-V','Thompson');