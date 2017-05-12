%% Create a Markov chain

N = 10; % number of states
A = 10; % number of actions

% create probabilities
Pr = rand(N,N,A);
Pr = Pr ./ repmat(sum(Pr,2),1,N);

% create rewards
R = rand(N,A); %rand(N,N,A);

% discount factor beta
beta = 0.99;

[V,pi] = value_iteration(Pr,R,beta);
[V2,pi2] = policy_iteration(Pr,R,beta);
assert(sum(abs(pi-pi2)) == 0);