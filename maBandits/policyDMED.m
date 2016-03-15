classdef policyDMED < Policy
    % The DMED policy of [Honda & Takemura, COLT 2010] in the special case of Bernoulli
    % rewards (can be used on any [0,1]-valued rewards, but warning: in the non-binary case, 
    % this is not the algorithm of [Honda & Takemura, COLT 2010])
    % (see note below on the variant)
    %
    % authors: Olivier Cappé, Aurélien Garivier

    % $Id: policyDMED.m,v 1.11 2012-06-05 13:26:38 cappe Exp $
        
    properties
        t % Index of the round
        lastAction % Stores the last action played
        N % Nb of times each action has been chosen
        S % Cumulated reward with each action
        L % Current list of actions to be played
        genuine = true;  % Variant: with this set to false, use a less aggressive list
                         % pruning criterion corresponding to the the version called
                         % DMED in [Garivier & Cappé, COLT 2011]; the default is the
                         % original proposal of [Honda & Takemura, COLT 2010] (called
                         % DMED+ in [Garivier & Cappé, COLT 2011])
    end
    
    methods
        function self = policyDMED(genuine)
            if (nargin == 1)
                self.genuine = genuine;
            end
        end
        function init(self, nbActions, horizon)
            self.t = 1;
            self.N = zeros(1, nbActions);
            self.S = zeros(1, nbActions);
            self.L = (1:nbActions);
        end
        function action = decision(self)
            if length(self.L) > 0
                % Play action on the list and update the list
                action = self.L(1);
                self.L = self.L(2:end);
            else
                % Make new list an play first action
                [~,c] = max(self.S./self.N); % Current best empirical mean
                if (self.genuine)
                    self.L = find(self.N.*kl(self.S./self.N, self.S(c)/self.N(c)) < ...
                      log(self.t./self.N));
                else
                    self.L = find(self.N.*kl(self.S./self.N, self.S(c)/self.N(c)) < ...
                      log(self.t));
                end
                action = self.L(1);
                self.L = self.L(2:end);
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

function y = kl(p,q)
  p = max(p, eps); p = min(p, 1-eps);  
  q = max(q, eps); q = min(q, 1-eps);
  y = p.*log(p./q) + (1-p).*log((1-p)./(1-q));
end
