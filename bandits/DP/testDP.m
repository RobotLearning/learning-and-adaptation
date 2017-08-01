%% Infinite-horizon MDP 
clc; clear; close all

N = 10; % number of states
A = 10; % number of actions

% create probabilities
Pr = rand(N,N,A);
Pr = Pr ./ repmat(sum(Pr,2),1,N);

% create rewards
R = rand(N,A); %rand(N,N,A);

% discount factor beta
beta = 0.99;

tic;
[V,pi] = value_iteration(Pr,R,beta);
t1 = toc
tic; 
[V2,pi2] = policy_iteration(Pr,R,beta);
t2 = toc
tic;
[V3,pi3] = lin_prog_mdp(Pr,R,beta);
t3 = toc
assert(sum(abs(pi-pi2)) == 0);
assert(sum(abs(pi2-pi3)) == 0);