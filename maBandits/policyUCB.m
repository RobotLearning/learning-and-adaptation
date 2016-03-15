classdef policyUCB < Policy
    % The classic UCB policy, with optimal tuning
    % of the constant (=1/2)
    %
    % Rewards are assumed to be bounded in [0,1]
    %
    % authors: Olivier Cappé and Aurélien Garivier
    % Ref.: [Auer, Cesa-Bianchi, Fischer, Machine Learning 2002; 
    %        Garivier & Cappé, COLT 2011]

    % $Id: policyUCB.m,v 1.5 2012-06-05 13:26:38 cappe Exp $

    properties
        c = 1/2 % Parameter of the UCB
        t % Number of the round
        lastAction % Stores the last action played
        N % Nb of times each action has been chosen
        S % Cumulated reward with each action
    end
    
    methods
        function self = policyUCB()
            
        end
        function init(self, nbActions, horizon)
            self.t = 1;
            self.N = zeros(1, nbActions);
            self.S = zeros(1, nbActions);
        end
        function action = decision(self)
            if any(self.N==0)
                action = find(self.N==0, 1);
            else
                ucb =  self.S./self.N + sqrt(self.c*log(self.t)./self.N);
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
