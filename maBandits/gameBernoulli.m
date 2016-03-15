classdef gameBernoulli<Game
    % Generic Bernoulli bandit game
    %
    % authors: Olivier Cappé and Aurélien Garivier

    % $Id: gameBernoulli.m,v 1.12 2012-06-05 13:26:38 cappe Exp $
    
    properties
        tabR % internal: array of all rewards
        N % internal: counters for rewards
    end
    
    methods
        function self = gameBernoulli(mu)
            self.mu = mu;
            self.nbActions = length(self.mu);
        end
        
        function [reward, action] = play(self, policy, n)
            policy.init(self.initRewards(n), n);
            reward = zeros(1, n);
            action = zeros(1, n);
            for t = 1:n
                action(t) = policy.decision();
                reward(t) = self.reward(action(t));
                policy.getReward(reward(t));
            end
        end
        
        function K = initRewards(self,n)
            % initiates the reward process, and returns the number of
            % actions
            K = length(self.mu);
            self.tabR = rand(K, n) < self.mu'*ones(1,n);
            self.N = zeros(1,K);            
        end
        
        function r = reward(self, a)
            self.N(a) = self.N(a) + 1;
            r = self.tabR(a, self.N(a));            
        end
    end    
end
