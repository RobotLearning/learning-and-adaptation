classdef gameExp<Game
    % (possibly truncated) exponential bandit game
    %
    % authors: Olivier Cappé and Aurélien Garivier

    % $Id: gameExp.m,v 1.8 2012-06-05 13:26:38 cappe Exp $

    properties
        param     % parameters of the Exponentials
        bound = 0 % if positive, the rewards are bounded exponential variables
                  %   that are normalized to [0,1]
        tabR      % internal: array of all rewards
        N         % internal: counters for rewards
    end
    
    methods
        function self = gameExp(param, bound)
            self.param = param;
            % compute expectations of the arms
            if (nargin >=2)
                self.bound = bound;
                self.mu = (1-exp(-self.param*bound))./self.param / bound;
            else 
                self.mu = 1./self.param;
            end
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
            %   actions
            K = length(self.param);
            self.tabR = -(1./self.param'*ones(1,n)) .* log(rand(K, n));
            % in the bounded case, rewards are normalized
            if self.bound>0
                self.tabR = min(self.tabR, self.bound)/self.bound;
            end
            self.N = zeros(1,K);            
        end
        
        function r = reward(self, a)
            self.N(a) = self.N(a) + 1;
            r = self.tabR(a, self.N(a));            
        end
    end
    
end
