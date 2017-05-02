%% Testing smooth regret

% Here we compare the regret and switching costs for Thompson
% and cautious thompson (or N-thompson)

clc; clear; close all;

% average over M experiments
M = 100;
% number of arms
K = 10;
% horizon, i.e. total num of time stages
N = 1000;

% scale for the switching cost
switch_alpha = 1/K;

num_algs = 2;
strategy{1}.name = 'Thompson-Cautious';
strategy{1}.a = 10;
strategy{1}.b = 1.0;
strategy{1}.sample_num = @(t) floor(t ^ (log(3*K/10)/log(N))); %floor(0.3*K);
strategy{2}.name = 'Thompson-Normal';
strategy{2}.a = 10;
strategy{2}.b = 1.0;

regret = zeros(num_algs,N);
cum_regret = zeros(num_algs,N);
switch_cost = zeros(num_algs,N);
cum_switch_cost = zeros(num_algs,N);
arms = zeros(num_algs,1);
last_arms = zeros(num_algs,1);

for j = 1:M % for each experiment
    
    % generating K reward parameters
    mu = 1*rand(K,1);
    var = 1*rand(K,1);
    
    % initialize different bandit strategies
    for k = 1:num_algs
        bandit{k} = Bandit(K,strategy{k});
    end
    
    for i = 1:N
        % generate rewards
        %rewards = mu + sqrt(var).*randn(K,1);
        % generate from uniform distribution
        rewards = (mu - sqrt(12*var)/2) + sqrt(12*var).*rand(K,1);
        
        % play bandit strategy
        for k = 1:num_algs
            arms(k) = bandit{k}.play();
            bandit{k}.reward(rewards(arms(k)),arms(k));
            % calculate regret
            regret(k,i) = max(mu) - rewards(arms(k));
            switch_cost(k,i) = switch_alpha * abs(arms(k) - last_arms(k));
            cum_regret(k,i:end) = cum_regret(k,i:end) + regret(k,i);
            cum_switch_cost(k,i:end) = cum_switch_cost(k,i:end) + switch_cost(k,i);
        end
        last_arms = arms;
        
    end
end

% plot cumulative regret
figure
subplot(2,1,1);
plot(cum_regret'/M);
legend('Cautious Thompson','Thompson','Location','northwest');
title('Growth of cumulative regret');
subplot(2,1,2);
plot((cum_regret + cum_switch_cost)'/M);
legend('Cautious Thompson','Thompson','Location','northwest');
title('Growth of cum switching cost + cum regret');
% plot switching costs
