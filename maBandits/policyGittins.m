classdef policyGittins < Policy
    % Gittins policy for a fixed horizon (in Bernoulli bandits)
    %
    % Note: To speed up computations this function uses a Matlab v7.3 mat
    % file of pre-computed values stored in matlab/GI.mat. Despite this file
    % being already big (30Mb) it contains only some data for horizons smaller
    % than 1000 and is not very dense for horizons larger than 500.        
    %
    % authors: Olivier Cappé, Aurélien Garivier, Emilie Kaufmann
    % Ref.: [Kaufmann, Cappé & Garivier, AISTATS 2012]

    % $Id: policyGittins.m,v 1.8 2012-06-06 09:42:04 cappe Exp $
    
    properties
        a,b % Parameters of the prior
        t % Index of the round
        lastAction % Stores the last action played
        N % Nb of times each action has been chosen
        S % Cumulated reward with each action
        horizon % Game length
        currentIndex % Current values of the Gittins indices
        GI % Pre-computed Gittins indices
    end
    
    methods
        function self = policyGittins(a,b)
            if nargin<2, b=1; end
            if nargin<1, a=1; end
            self.a = a;
            self.b = b;
            load('data/GI.mat');
            self.GI = GI;
        end
        
        function init(self, nbActions, horizon)
            self.t = 1;
            self.N = zeros(1, nbActions);
            self.S = zeros(1, nbActions);
            self.horizon = horizon;
            if (self.horizon > length(self.GI))
                warning('Horizon larger than %d, slow computations ahead',  length(self.GI));
            end
            self.currentIndex = ones(1, nbActions);
            self.lastAction = 1;
        end
        
        function action = decision(self)
           % Make decision based on Gittins indices
           K = length(self.N);
           newIndex = self.gittinsIndex(self.a+self.S(self.lastAction), self.b+self.N(self.lastAction)-self.S(self.lastAction), self.horizon-self.t+1);
           % If current index is increased, then play the same arm and avoid
           % computation of the other indices
           if newIndex > self.currentIndex(self.lastAction)
               action = self.lastAction;
               self.currentIndex(action) = newIndex; % Other indices are not updated but can
                                                     % only decrease
           else
               self.currentIndex(self.lastAction) = newIndex;
               for i = [1:self.lastAction-1, self.lastAction+1:K]
                   self.currentIndex(i) = self.gittinsIndex(self.a+self.S(i), self.b+self.N(i)-self.S(i), self.horizon-self.t+1);
               end
               m = max(self.currentIndex); I = find(self.currentIndex == m);
               action = I(1+floor(length(I)*rand));
           end
           self.lastAction = action;
        end

        function x = gittinsIndex(self, a, b, T)
            % Actual Gittins index computation (if time to horizon is
            % small, will try to use precomputed values)
            if T<=length(self.GI) && a<=size(self.GI{T}, 1) && b<=size(self.GI{T}, 2) && self.GI{T}(a,b)>0, x = self.GI{T}(a,b);
            else
                % Faster search, using convexity
                tol = 1e-7;
                l = a/(a+b);
                c = policyGittins.gameValue(l,a,b,T);
                if c<=tol, x = l;
                else
                    u = 1;
                    d = policyGittins.gameValue(u,a,b,T);
                    while(d==0), u = (l+u)/2; d = policyGittins.gameValue(u,a,b,T); end
                    x=u;
                    while d>tol
                        nx = x + d * (x-l) / (c-d); 
                        l = x; c = d; 
                        x = nx; d = policyGittins.gameValue(x,a,b,T); % >0 theoretically, =0 if value was linear on [l,x]
                    end
                end
                self.GI{T}(a,b) = x;
            end
        end

        function getReward(self, reward)
            self.N(self.lastAction) = self.N(self.lastAction) + 1; 
            self.S(self.lastAction) = self.S(self.lastAction)  + reward;
            self.t = self.t + 1;
        end
        
        function saveGI(self)
            % To save the pre-computed table (add to calls to this function is you
            % wantr to update your own table of pre-computed Gittins indices
            GI = self.GI;
            save('data/GI.mat', 'GI');
        end
        
        function action = decision2(self)
           % Alternative version: indices are all re-computed at each step
           K = length(self.N);
           self.currentIndex = zeros(1, K);
           for i = 1:K
               self.currentIndex(i) = self.gittinsIndex(self.a+self.S(i), self.b+self.N(i)-self.S(i), self.horizon-self.t+1);
           end
           m = max(self.currentIndex); I = find(self.currentIndex == m);
           action = I(1+floor(length(I)*rand));
           self.lastAction = action;
        end
                
        function x = gittinsIndex2(self, a, b, T)
            % Alternative version: slower search using dichotomy, safer but slower
            if T<=length(self.GI) && a<=size(self.GI{T}, 1) && b<=size(self.GI{T}, 2) && self.GI{T}(a,b)>0, x = self.GI{T}(a,b);
            else
                l = a/(a+b);
                u = 1;
                while(u-l>1e-5)
                    x = (l+u)/2;
                    if policyGittins.gameValue(x,a,b,T)==0, u = x; else l=x; end
                end
                self.GI{T}(a,b) = x;
            end
        end
                 
    end
    
    methods(Static)
        % Internal functions used for computing Gittins indices
        function v = gameValue(l,a,b,T)
            V = policyGittins.snellEnv(l,a,b,T);
            v = V(1,1);
        end

        function [V, A] = snellEnv(l,a,b,T)
            V = zeros(T); % V(t,j) = value of the game at time i, given that sum(X(1:(t-1)))=j-1
            A = zeros(T); % A(t,j) = value of the previsible part of the snell enveloppe at time t+1, given S(t)=j-1 
            V(:, T) = (1:T)-1 - (T-1)*l;
            for s=(T-1):-1:1
                for j = 1:s
                    p = (j-1+a) / (s-1+a+b);
                    V(j, s) = max(j-1-(s-1)*l, p*V(j+1, s+1) + (1-p)*V(j, s+1));
                    A(j, s) = max(j-1-(s-1)*l - (p*V(j+1, s+1) + (1-p)*V(j, s+1)), 0);
                end
            end
        end        
    end

end
