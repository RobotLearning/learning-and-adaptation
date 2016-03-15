classdef gamePoisson<Game
    % (possibly truncated) Poisson bandit game
    %
    % authors: Olivier Cappé and Aurélien Garivier

    % $Id: gamePoisson.m,v 1.5 2012-06-05 13:26:38 cappe Exp $
    
    properties
        param     % parameters of the Poisson
        bound = 0 % if positive, the rewards are bounded poisson variables
                  %   that are normalized to [0,1]
        tabR      % internal: array of all rewards
        N         % internal: counters for rewards
    end
    
    methods
        function self = gamePoisson(param, bound)
            self.param = param;
            self.nbActions = length(self.param);
            if (nargin>=2)
                self.bound = bound;
                % Compute expectations of the arms
                self.mu = zeros(1, self.nbActions);
                for a = 1:self.nbActions
                    p = exp(-param(a)) * param(a).^(0:(bound-1)) ./ factorial(0:(bound-1));
                    p = [p, 1-sum(p)];
                    self.mu(a) = p * (0:bound)' / bound;
                end
            else
                self.mu = self.param;
            end
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
            K = length(self.param);
            for a=1:K
                self.tabR(a, :) = poissrnd(self.param(a), 1, n);
            end
            % in the bounded case, rewards are normalized
            if self.bound>0
                self.tabR = min(self.tabR, self.bound) / self.bound;
            end
            self.N = zeros(1,K);            
        end
        
        function r = reward(self, a)
            self.N(a) = self.N(a) + 1;
            r = self.tabR(a, self.N(a));            
        end
    end
    
end
