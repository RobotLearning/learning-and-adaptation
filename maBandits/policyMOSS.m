classdef policyMOSS < Policy
    % MOSS policy
    %
    % Rewards are assumed to be bounded in [0,1]
    %
    % authors: Olivier Cappé and Aurélien Garivier
    % Ref.: [Audibert & Bubeck, JMLR 2010]

    % $Id: policyMOSS.m,v 1.5 2012-06-05 13:26:38 cappe Exp $
    
    properties
        t % number of the round
        horizon
        lastAction % stores the last action played
        N % nb of times each action has been chosen
        S % cumulated reward with each action
    end
    
    methods
        function self = policyMOSS()
            
        end
        function init(self, nbActions, horizon)
            self.horizon = horizon; % MOSS en a besoin
            self.t = 1;
            self.N = zeros(1, nbActions);
            self.S = zeros(1, nbActions);
        end
        function action = decision(self)
            if any(self.N==0)
                action = find(self.N==0, 1);
            else
                ucb =  self.S./self.N + sqrt(log(max(self.horizon./(length(self.N)*self.N)))./self.N);
                m = max(ucb); I = find(ucb == m);
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
