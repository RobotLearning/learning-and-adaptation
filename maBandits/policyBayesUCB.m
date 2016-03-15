classdef policyBayesUCB < Policy
    % Bayes-UCB for Bernoulli bandits
    % (using uniform prior and 1-1/t quantile)
    % 
    % authors: Olivier Cappé, Aurélien Garivier
    % Ref.: [Kaufmann, Cappé & Garivier, AISTATS 2012]

    % $Id: policyBayesUCB.m,v 1.6 2012-06-05 13:26:38 cappe Exp $

    properties
        t % Number of the round
        lastAction % Stores the last action played
        N % Number of times each action has been chosen
        S % Cumulated reward with each action

    end
    
    methods
        function self = policyBayesUCB()
            
        end
        
        function init(self, nbActions,horizon)
            self.t = 1;
            self.N = zeros(1, nbActions);
            self.S = zeros(1, nbActions);
        end
        
        function [action,quantiles] = decision(self)
            if any(self.N==0)
                action = find(self.N==0, 1);
                quantiles = ones(size(self.N));
            else
                quantiles = betaincinv(1-1/(self.t),self.S+1,self.N-self.S + 1);
                % Using beta quantile
                m = max(quantiles); I = find(quantiles == m);
                action = I(1+floor(length(I)*rand));
            end
            self.lastAction = action;
        end
        
        function getReward(self, reward)
            self.N(self.lastAction) = self.N(self.lastAction) + 1; 
            self.S(self.lastAction) = self.S(self.lastAction)  + reward;
            self.t = self.t + 1;
        end        
    end
    
end