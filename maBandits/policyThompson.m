classdef policyThompson < Policy
    % Thompson policy for Bernoulli bandits
    %
    % authors: Olivier Cappé, Aurélien Garivier and Emilie Kaufmann
    % Ref.: [Kaufmann, Korda & Munos, arXiv:1205.4217]

    % $Id: policyThompson.m,v 1.5 2012-06-05 13:26:38 cappe Exp $
    
    properties
        t % Number of round
        lastAction % Stores the last action played
        N % Nb of times each action has been chosen
        S % Cumulated reward with each action
    end
    
    methods
        function self = policyThompson()
            
        end
        function init(self, nbActions, horizon)
            self.t = 1;
            self.N = zeros(1, nbActions);
            self.S = zeros(1, nbActions);
        end
        function action = decision(self)
            % Draws from posteriors
            u = betarnd(1+self.S, 1+self.N-self.S);
            [m, action] = max(u);
            self.lastAction = action;
        end
        function getReward(self, reward)
            self.N(self.lastAction) = self.N(self.lastAction) + 1; 
            self.S(self.lastAction) = self.S(self.lastAction) + reward;
            self.t = self.t + 1;
        end        
    end

end
