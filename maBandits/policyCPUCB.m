classdef policyCPUCB < Policy
    % Copper-Pearson UCB 
    %
    % authors: Olivier Cappé, Aurélien Garivier
    % Ref.: [Garivier & Cappé, COLT 2011]

    % $Id: policyCPUCB.m,v 1.6 2012-06-05 13:26:38 cappe Exp $
        
    properties
        c = 1.01 % parameter of the ucb
        t % number of the round
        lastAction % stores the last action played
        N % nb of times each action has been chosen
        S % cumulated reward with each action
    end
    
    methods
        function self = policyCPUCB()
            
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
                K = length(self.N); 
                [~, ic] = binofit(self.S, self.N, 1/(self.t^self.c));
                ucb = ic(:,2);
                m = max(ucb); I = find(ucb == m);
                action = I(1+floor(length(I)*rand));
            end
            self.lastAction = action;
        end
        function getReward(self, reward)
            self.N(self.lastAction) = self.N(self.lastAction) + 1; 
            self.S(self.lastAction) = self.S(self.lastAction) + reward;
            self.t = self.t + 1;
        end
    end
    
end
