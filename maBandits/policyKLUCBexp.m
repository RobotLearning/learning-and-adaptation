classdef policyKLUCBexp < Policy
    % The KL-UCB policy for Exponential rewards
    %
    % authors: Olivier Cappé and Aurélien Garivier
    % Ref.: [Garivier & Cappé, COLT 2011]

    % $Id: policyKLUCBexp.m,v 1.11 2012-06-05 13:26:38 cappe Exp $
    
    properties
        c = 0 % parameter of the ucb
        t % number of the round
        lastAction % stores the last action played
        N % nb of times each action has been chosen
        S % cumulated reward with each action
        variant = false % Variant: with this set to true, use a slightly
                         % lower heuristic UCB (that is, the version called
                         % KLUCB+ in [Garivier & Cappé, COLT 2011])
       bound = 1 % if bounded is not 1, then it is assumed that 
                   % the rewards are actually bounded expontential
                   % variables that have been normalized to [0,1] (i.e.,
                   % divided by bound)
    end
    
    methods
        function self = policyKLUCBexp(c, variant, bound)
            if (nargin >= 1)
                self.c = c;
            end
            if (nargin >= 2)
                self.variant = variant;
            end     
            if (nargin >= 3)
                self.bound = bound;
            end                        
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
                if (self.variant)
                    ucb = klIC(self.S./self.N, ...
                      (log(self.t./self.N) ...
                      + self.c*log(log(self.t)))./self.N);
                else
                    ucb = klIC(self.S./self.N, ...
                      (log(self.t) + self.c*log(log(self.t)))./self.N);
                end
                m = max(ucb); I = find(ucb == m);
                action = I(1+floor(length(I)*rand));
            end
            self.lastAction = action;
        end
        
        function getReward(self, reward)
            self.N(self.lastAction) = self.N(self.lastAction) + 1; 
            self.S(self.lastAction) = self.S(self.lastAction) + self.bound * reward;
            self.t = self.t + 1;
        end
    end
end

% Vectorised search of the UCB of all arms by dichotomy
function q = klIC(p, d)
    l = p; 
    petit = d<1/3;
    u(petit) = p(petit) ./( 1 - sqrt(2*d(petit)));
    u(~petit)= p(~petit) .* exp(1+d(~petit));
    while max(u-l)>1e-4
        m = (u+l)/2;
        down = kl(p,m) > d;
        u(down) = m(down);
        l(~down) = m(~down);
    end
    q = u;
end

function y = kl(p,q)
  p = max(p, eps);
  q = max(q, eps);
  y = p./q - 1 - log(p./q); 
end

