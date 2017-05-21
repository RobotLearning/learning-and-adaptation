%% Testing smooth regret

% Here we compare the regret and switching costs for Thompson
% and cautious thompson (or N-thompson)

clc; clear; close all;

rng(13);
% average over M experiments
M = 1; %100;
% number of arms
K = 20;
% horizon, i.e. total num of time stages
N = 500;

% scale for the switching cost
switch_alpha = 1/K;

num_algs = 3;
strategy{1}.name = 'Thompson-Cautious';
strategy{1}.a = 10;
strategy{1}.b = 2.0;
strategy{1}.sample_num = @(t) floor(0.3*K); %0.3*K 
strategy{2}.name = 'Thompson-Cautious';
strategy{2}.a = 10;
strategy{2}.b = 2.0;
strategy{2}.sample_num = @(t) floor(0.5*K); %0.3*K 
strategy{3}.name = 'Thompson-Normal';
strategy{3}.a = 10;
strategy{3}.b = 2.0;
% strategy{3}.name = 'Thompson-Regularized';
% strategy{3}.a = 10;
% strategy{3}.b = 2.0;
% strategy{3}.lambda = switch_alpha/20;
% strategy{4}.name = 'EPS-GREEDY';
% strategy{4}.eps = @(t) 10/t;
% strategy{5}.name = 'UCB';
% strategy{5}.rho = 1;
% strategy{6}.name = 'UCBV';
% strategy{7}.name = 'UCBV-Regularized'; % regularized UCBV
% strategy{7}.lambda = switch_alpha/20;

regret = zeros(num_algs,N);
cum_regret = zeros(num_algs,N);
switch_cost = zeros(num_algs,N);
cum_switch_cost = zeros(num_algs,N);
arms = zeros(num_algs,N);
last_arms = zeros(num_algs,1);

for j = 1:M % for each experiment
    
    fprintf('Experiment %d\n', j);
    
    % generating K reward parameters
    mu = 1*rand(K,1);
    var = 1*rand(K,1);
    
    % initialize different bandit strategies
    for k = 1:num_algs
        band{k} = Bandit(K,strategy{k});
    end
    
    for i = 1:N
        % generate rewards
        rewards = mu + sqrt(var).*randn(K,1);
        % generate from uniform distribution
        %rewards = (mu - sqrt(12*var)/2) + sqrt(12*var).*rand(K,1);
        
        % play bandit strategy
        for k = 1:num_algs
            arms(k,i) = band{k}.play();
            band{k}.reward(rewards(arms(k,i)),arms(k,i));
            % calculate regret
            regret(k,i) = max(mu) - rewards(arms(k,i));
            switch_cost(k,i) = switch_alpha * abs(arms(k,i) - last_arms(k));
            cum_regret(k,i:end) = cum_regret(k,i:end) + regret(k,i);
            cum_switch_cost(k,i:end) = cum_switch_cost(k,i:end) + switch_cost(k,i);
        end
        last_arms = arms(:,i);
        
    end
end

% % plot cumulative regret
% figure
% subplot(2,1,1);
% plot(cum_regret'/M);
% legend(strategy{1}.name,...
%     strategy{2}.name,...
%     strategy{3}.name,...
%     strategy{4}.name,...
%     strategy{5}.name,...
%     strategy{6}.name,...
%     strategy{7}.name,...
%        'Location','northwest');
% title('Average cumulative regret');
% subplot(2,1,2);
% plot((cum_regret + cum_switch_cost)'/M);
% legend(strategy{1}.name,...
%     strategy{2}.name,...
%     strategy{3}.name,...
%     strategy{4}.name,...
%     strategy{5}.name,...
%     strategy{6}.name,...
%     strategy{7}.name,...
%        'Location','northwest');
% title('Average cumulative augmented regret');
% % plot switching costs

figure;
plot(arms(3,1:50),'-b','LineWidth',1.0);
hold on
plot(arms(1,1:50),'-r','LineWidth',1.0);
plot(arms(2,1:50),'--r','LineWidth',1.0);
legend('Thompson-Normal',...
       'Thompson Cautious (S = 3)', ...
       'Thompson Cautious (S = 5)');
title('Switching between arms');
axis tight