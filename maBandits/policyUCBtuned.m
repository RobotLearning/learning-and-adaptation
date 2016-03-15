classdef policyUCBtuned < Policy
    % UCB tuned (with variance estimation)
    %
    % authors: Olivier Cappé and Aurélien Garivier
    % Ref.: [Auer, Cesa-Bianchi, Fischer, Machine Learning 2002]

    % $Id: policyUCBtuned.m,v 1.7 2012-06-05 13:32:34 cappe Exp $
    
    properties
        t % Index of the round
        lastAction % Stores the last action played
        N % Nb of times each action has been chosen
        S % Cumulated reward with each action
        S2 % Cumulated squared reward with each action
    end
    
    methods
        function self = policyUCBtuned()
            
        end
        function init(self, nbActions, horizon)
            self.t = 1;
            self.N = zeros(1, nbActions);
            self.S = zeros(1, nbActions);
            self.S2 = zeros(1, nbActions);
        end
        function action = decision(self)
            if any(self.N==0)
                action = find(self.N==0, 1);
            else
                m =  self.S./self.N;
                V = self.S2./self.N - m.^2 + sqrt(2*log(self.t)./self.N); % Note: this
                                                                          % correction
                                                                          % makes sense
                                                                          % for rewards
                                                                          % in [0,1]
                %ucb = m + sqrt(log(self.t)./self.N.*min(1/4, V)) % This one could also be used
                ucb = m + sqrt(log(self.t)./self.N.*V);
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
