classdef policyKLUCB < Policy
    % The KL-UCB policy for bounded variables
    % (based on binary Kullback-Leibler divergence)
    %
    % authors: Olivier Cappé and Aurélien Garivier
    % Ref.: [Garivier & Cappé, COLT 2011]

    % $Id: policyKLUCB.m,v 1.20 2012-06-05 13:26:38 cappe Exp $
        
    properties
        c = 0 % parameter of the ucb
        t % number of the round
        lastAction % stores the last action played
        N % nb of times each action has been chosen
        S % cumulated reward with each action
        variant = false % Variant: with this set to true, use a slightly
                         % lower heuristic UCB (that is, the version called
                         % KLUCB+ in [Garivier & Cappé, COLT 2011])
    end
    
    methods
        function self = policyKLUCB(c, variant)
            if (nargin >= 1)
                self.c = c;
            end
            if (nargin >= 2)
                self.variant = variant;
            end
        end
        
        function init(self, nbActions, horizon)
            self.t = 1;
            self.N = zeros(1, nbActions);
            self.S = zeros(1, nbActions);
        end        
        
        function [action, ucb] = decision(self)
            if any(self.N==0)
                action = find(self.N==0, 1);
                ucb = ones(size(self.N));
            else
                if (self.variant)
                    ucb = klIC(self.S./self.N, (log(self.t./self.N) ...
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
            self.S(self.lastAction) = self.S(self.lastAction) + reward;
            self.t = self.t + 1;
        end
    end
end

% Vectorised search of the UCB of all arms by dichotomy
function qM = klIC(p, d)
    lM = p; uM = min(1,p+sqrt(d/2)); % ones(size(p));
    for j = 1:16
        qM = (uM+lM)/2;
        down = kl(p,qM) > d;
        uM(down) = qM(down);
        lM(~down) = qM(~down);
    end
    qM = uM;
end

function y = kl(p,q)
  p = max(p, eps); p = min(p, 1-eps);  
  q = max(q, eps); q = min(q, 1-eps);
  y = p.*log(p./q) + (1-p).*log((1-p)./(1-q));
end

% This is much longer (using fzero)...
function qM = klIC_alt(p, d)
    if (kl_shift(p,1,d) <= 0)
        qM = 1;
    else
        lM = p; uM = min(1,p+sqrt(d/2));
        options = optimset('TolX',1e-5);
        qM = fzero(@(q) kl_shift(p,q,d), [lM uM], options);
    end
end

function y = kl_shift(p,q,d)
  y = kl(p,q)-d;
end
