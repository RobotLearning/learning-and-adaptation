classdef policyKLempUCB < Policy
    % The KL-UCB policy, based on Empirical Likelihood upper-confidence
    % bounds
    %
    % Rewards are assumed to be bounded in [0,1]
    %
    % authors: Olivier Cappé and Aurélien Garivier
    % Ref.: [Cappé, Garivier, Maillard, Munos, Stoltz, 2012]

    % $Id: policyKLempUCB.m,v 1.8 2012-06-05 13:28:37 cappe Exp $
    
    properties
        c = 0 % parameter of the ucb
        t % number of the round
        lastAction % stores the last action played
        X % X{a} = vector of rewards seen on arm a
        Nb % Nb{a}(x) = number of times reward X{a}(x) has been seen on arm a
        variant = false % Variant: with this set to true, use a slightly
                        % lower heuristic UCB
    end
    
    methods
        function self = policyKLempUCB(c, variant)
            if (nargin >= 1)
                self.c = c;
            end
            if (nargin >= 2)
                self.variant = variant;
            end            
        end
        function init(self, nbActions, horizon)
            self.t = 1;
            self.X = cell(1, nbActions);
            self.Nb = cell(1, nbActions);
        end        
        function [action, ucb] = decision(self)
            K = length(self.Nb); 
            ucb = ones(1, K);
            if any(cellfun(@isempty, self.X))
                action = find(cellfun(@isempty, self.X), 1);
            else
                for a=1:K
                    if (self.variant)
                        ucb(a) = KLucb(self.X{a}, self.Nb{a}, (log(self.t/self.Nb{a}) + self.c*log(log(self.t)))./sum(self.Nb{a}));
                    else
                        ucb(a) = KLucb(self.X{a}, self.Nb{a}, (log(self.t) + self.c*log(log(self.t)))./sum(self.Nb{a}));
                    end
                end
                m = max(ucb); I = find(ucb == m);
                action = I(1+floor(length(I)*rand));
            end
            self.lastAction = action;
        end
        function getReward(self, reward)
            k = find(self.X{self.lastAction}==reward);
            if isempty(k)
                self.X{self.lastAction} = [self.X{self.lastAction}, reward]; 
                self.Nb{self.lastAction} = [self.Nb{self.lastAction}, 1]; 
            else
                self.Nb{self.lastAction}(k) = self.Nb{self.lastAction}(k) + 1;
            end
            self.t = self.t + 1;
        end
        function n = N(self)
            n = cellfun(@sum,self.Nb);
        end        
    end
end


function [u,q] = KLucb(x, n, d)
    v = x';
    p = n / sum(n);
    if all(x ~= 1)
        v = [v; 1];
        p = [p, 0];
    end
    q = maxEV(p, v, d);
    u = q*v;
end


function Uq = maxEV(p, V, klMax)
    % Maximize expectation of V wrt. q st. KL(p,q) < klMax.
    % Ref.:  Section 3.2 of [Filippi, Cappé & Garivier - Allerton, 2011].
    
    k = length(p); Uq = zeros(1, k);
    K = find(p==0);
    Kb = find(p>0); 
    if ~isempty(K) % some components of p are zero
        J = find(V(K) == max(V(K))); % argmax_{i\in K} V(i)
        eta = V(K(J(1))); % max_{i\in K} V(i) 
        y = p(Kb)*log(eta-V(Kb)) + log(p(Kb) * (1./(eta-V(Kb))));
        if (eta > max(V(Kb))) && ( y < klMax)
            r = 1-exp( y-klMax);
            Uqtemp = p(Kb)./(eta - V(Kb))';
            Uq(Kb) = Uqtemp*(1-r)/sum(Uqtemp);
            Uq(K(J(1))) = r; % or Uq(K(J)) = r / length(J);
            return;
        end
    end
    % Here, only points where p is strictly positive (in Kb) will receive
    % non-zero mass.
    if any(abs(V(Kb) - V(Kb(1)))>1e-8)
        eta = reseqp(p(Kb), V(Kb), klMax); 
        Uq = p./(eta-V)';
        Uq=Uq/sum(Uq);
    else
        Uq(Kb) = 1/length(Kb); % Case where all values in V(Kb) are almost identical.
    end
end


function l = reseqp(p, V, klMax)
    % Solve f(reseqp(p, V, klMax)) = klMax using Newton's method.
    % Note: This is a subroutine of maxEV.
    % Reference: Eq. (4) in Section 3.2 of [Filippi, Cappé & Garivier - Allerton, 2011].
    
    mV = max(V);
    l = mV+0.1; tol = 1e-4;
    if mV<min(V)+tol, l = inf; return; end
    u = p * (1./(l-V));
    y = p*log(l-V) + log(u) - klMax; 
    while abs(y)>tol
        yp = u -  p * (1./(l-V).^2) / u; % derivative
        l = l - y / yp; % newton iteration
        if l<mV, l = (l+y/yp +mV)/2; end % unlikely, but not impossible
        u = p * (1./(l-V));
        y = p*log(l-V) + log(u) - klMax;  
    end
end
