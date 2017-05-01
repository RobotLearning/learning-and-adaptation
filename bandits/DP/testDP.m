%% Create a Markov chain

N = 2; % number of states
A = 2; % number of actions

% create probabilities
Pr = rand(N,N,A);
Pr = Pr ./ sum(Pr,2);

% create rewards
R = rand(N,N,A);

% discount factor beta
beta = 0.99;

[V,pi] = value_iteration(Pr,R,beta);
[V2,pi2] = policy_iteration(Pr,R,beta);