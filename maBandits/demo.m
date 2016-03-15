% Script to illustrate the various functions in maBandits
%
% authors: Olivier Cappé, Aurélien Garivier, Emilie Kaufmann

% $Id: demo.m,v 1.33 2012-06-05 13:26:38 cappe Exp $


%%% First scenario: Bernoulli bandits
disp('--- First scenario: Bernoulli bandits');

% Scenario of Figure 1 in [Garivier & Cappé, COLT 2011]
% mu = [0.9, 0.8]; game = gameBernoulli(mu); fname = 'results/2a';

% Scenario of Figure 2 in [Garivier & Cappé, COLT 2011]
% mu = [0.1 0.05 0.05 0.05 0.02 0.02 0.02 0.01 0.01 0.01]; game = gameBernoulli(mu); fname = 'results/10';

% Scenario of Figure 1 in [Kaufmann, Cappé & Garivier, AISTATS 2012]
mu = [0.2, 0.1]; game = gameBernoulli(mu); fname = 'results/2abis';
% Choice of policies to be run
policies = {policyGittins(), policyBayesUCB(), policyThompson(), policyKLUCB(), ...
            policyCPUCB(), policyUCB(), policyUCBtuned(), policyUCBV(), policyMOSS(), ...
            policyDMED(), policyKLempUCB};
% n is length of play, N is number of plays 
n = 500; N = 20;
% Time indices at which the results will be saved
tsave = round(linspace(1,n,200));

% Run everything one policy after each other
defaultStream = RandStream.getDefaultStream; 
savedState = defaultStream.State;
for k = 1:length(policies)
    defaultStream.State = savedState;
    tic; experiment(game, n, N, policies{k}, tsave, fname); toc 
end
plotResults(game, n, N, policies, fname);


disp('Type any key to continue...'); pause;         
       

%%% Second scenario: Bounded Exponential rewards
disp('--- Second scenario: Bounded Exponential rewards');

% Scenario of Figure 3 in [Garivier & Cappé, COLT 2011]
param = 1./[5 4 3 2 1]; game = gameExp(param, 10); fname = 'results/E5b';
policies = {policyKLUCB(), policyUCB(), policyUCBtuned(), policyUCBV(), policyMOSS(), ...
           policyKLempUCB, policyKLUCBexp(0, false, 10)};
n = 2000; N = 20;
tsave = round(linspace(1,n,200));
         
defaultStream = RandStream.getDefaultStream; 
savedState = defaultStream.State;
for k = 1:length(policies)
    defaultStream.State = savedState;
    tic; experiment(game, n, N, policies{k}, tsave, fname); toc 
end
plotResults(game, n, N, policies, fname);


disp('Type any key to continue...'); pause;         


%%% Third scenario: (unbounded) Exponential rewards
disp('--- Third scenario: (unbounded) Exponential rewards');

% Unbounded version of the scenario of Figure 3 in [Garivier & Cappé, COLT 2011]
param = 1./[5 4 3 2 1]; game = gameExp(param); fname = 'results/E5';
policies = {policyUCBtuned(), policyKLUCBexp()}; % warning: the parameters of
                                                 % policyUCBtuned, tuned for
                                                 % [0,1]-valued rewards, 
                                                 % are left unchanged in this experiment
n = 2000; N = 20;
tsave = round(linspace(1,n,200));
         
defaultStream = RandStream.getDefaultStream; 
savedState = defaultStream.State;
for k = 1:length(policies)
    defaultStream.State = savedState;
    tic; experiment(game, n, N, policies{k}, tsave, fname); toc 
end
plotResults(game, n, N, policies, fname);


disp('Type any key to continue...'); pause;         


%%% Fourth scenario: Bounded Poisson rewards
disp('--- Fourth scenario: Bounded Poisson rewards');

% Scenario of Figure 3 in [Cappé, Garivier, Maillard, Munos, Stoltz, 2012]
param = 0.5+(1:6)/3; game = gamePoisson(param, 10); fname = 'results/P6b';
policies = {policyKLUCB(), policyUCB(), policyUCBtuned(), policyUCBV(), policyMOSS(), ...
            policyKLempUCB, policyKLUCBpoisson(0, false, 10)};
n = 2000; N = 20;
tsave = round(linspace(1,n,200));
         
defaultStream = RandStream.getDefaultStream; 
savedState = defaultStream.State;
for k = 1:length(policies)
    defaultStream.State = savedState;
    tic; experiment(game, n, N, policies{k}, tsave, fname); toc 
end
plotResults(game, n, N, policies, fname);
