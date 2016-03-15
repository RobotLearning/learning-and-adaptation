classdef policyUCBV < Policy
    % UCBV
    %
    % Rewards are assumed to be bounded in [0,1]
    %
    % authors: Olivier Cappé and Aurélien Garivier
    % Ref.: [Audibert, Munos & Szepesvari, 2009]

    % $Id: policyUCBV.m,v 1.8 2012-06-05 13:26:38 cappe Exp $
    
    properties
        t % index of the round
        lastAction % stores the last action played
        N % nb of times each action has been chosen
        S % cumulated reward with each action
        S2 % cumulated squared reward with each action
    end
    
    methods
        function self = policyUCBV()
            
        end
        function init(self, nbActions, horizon)
            self.t = 1;
            self.N = zeros(1, nbActions);
            self.S = zeros(1, nbActions);
            self.S2 = zeros(1, nbActions);
        end
        function action = decision(self)
            if any(self.N<=1)
                action = find(self.N<=1, 1);
            else
                m =  self.S./self.N;
                V = self.S2./self.N - m.^2;
                ucb = m + sqrt(2*log(self.t).*V./self.N) + 3*log(self.t)./self.N;
                m = max(ucb); I = find(ucb == m);
                action = I(1+floor(length(I)*rand));
            end
            self.lastAction = action;
        end
        function getReward(self, reward)
            self.N(self.lastAction) = self.N(self.lastAction) + 1; 
            self.S(self.lastAction) = self.S(self.lastAction)  + reward;
            self.S2(self.lastAction) = self.S2(self.lastAction)  + reward^2;
            self.t = self.t + 1;
        end        
    end
    
end
