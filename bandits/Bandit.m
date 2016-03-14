%% Bandit class for implementing different strategies

classdef Bandit < handle
    
    properties
       
        % strategies for exploring/exploiting
        strategy
        % mean and variance statistics for each arm
        mean
        var
        % counters for number of plays so far
        cnt
        
    end
    
    methods
        
        %% Initialize Bandit with a strategy
        function obj = Bandit(arms,str)
           
            obj.mean = zeros(arms,1);
            obj.var = zeros(arms,1);
            obj.strategy = str;
            obj.cnt = zeros(arms,1);
            
            if strcmp(obj.strategy.name,'UCB1-N') || ...
                strcmp(obj.strategy.name,'UCB1-V')
                obj.strategy.q = zeros(arms,1); % sum of sqr rewards
            end
            
        end
        
        %% Strategies for playing
        function idx = play(obj)
            
            arms = length(obj.mean);
            switch obj.strategy.name
                case 'EPS-GREEDY'
                    idx = obj.epsGreedy(arms);                    
                case 'EPS-FIRST'
                    idx = obj.epsFirst(arms);
                case 'UCB1'
                    idx = obj.ucb1(arms);
                case 'UCB1-N'
                    idx = obj.ucb1Normal(arms);
                case 'UCB1-V'
                    idx = obj.ucb1Variance(arms);
                otherwise
                    error('Alg not implemented!');
            end           
            
        end
        
        % classical ucb strategy
        function idx = ucb1(obj,arms)
            
            ntot = sum(obj.cnt) + 1;
            % make sure each machine is played once
            rho = obj.strategy.rho;
            % for variance aware ucb
%             qs = obj.strategy.q ./ obj.var;
%             V = qs - obj.mean.^2 + sqrt(2*log(obj.cnt)./obj.var);
%             rho = min(0.25*ones(length(V),1),V);
            if ntot <= arms
                idx = ntot;
            else
                [~,idx] = max(obj.mean + sqrt(rho.*log(ntot))./sqrt(obj.cnt));
            end
                
        end
        
        % classical ucb with variance estimate
        function idx = ucb1Variance(obj,arms)
            
            ntot = sum(obj.cnt) + 1;
            % for variance aware ucb
            qs = obj.strategy.q ./ obj.var;
            V = qs - obj.mean.^2 + sqrt(2*log(ntot)./obj.cnt);
            rho = min(0.25*ones(length(V),1),V);
            if ntot <= arms
                idx = ntot;
            else
                [~,idx] = max(obj.mean + sqrt(rho.*log(ntot))./sqrt(obj.cnt));
            end
        end
        
        % classical ucb assuming normal distribution
        function idx = ucb1Normal(obj,arms)
            
            ntot = sum(obj.cnt) + 1;
            n = obj.cnt;
            val = 1; % not 8
            if any(n < ceil(val*log(ntot)))
                idx = find(n < ceil(val*log(ntot)));
                idx = idx(randi(length(idx)));
            else
                qs = obj.strategy.q;
                val = log(ntot-1) * (qs - n.*obj.mean.^2) ./ (n.*(n-1));
                [~,idx] = max(obj.mean + 4*sqrt(val));
            end
        end
        
        % eps-n greedy with n the number of time steps
        % decaying exploration rate
        function idx = epsGreedy(obj,arms)
            
            ntot = sum(obj.cnt) + 1;
            eps = rand();
            maxval = obj.strategy.eps(ntot);
            if eps < maxval %obj.strategy.eps 
                % explore randomly
                idx = randi(arms);
            else
                % greedy action
                [~,idx] = max(obj.mean);
            end 
            
        end
        
        function idx = epsFirst(obj,arms)
            
            ntot = sum(obj.cnt) + 1;
            if ntot < obj.strategy.eps * obj.strategy.N 
                % sample at random
                idx = randi(arms);
            else
                % greedy action
                [~,idx] = max(obj.mean);
            end
        end
        
        %% Get rewards and update (sufficient) statistics
        function reward(obj,R,I)
            
            obj.mean(I) = (obj.mean(I)*obj.cnt(I) + R)/(obj.cnt(I)+1);
            obj.cnt(I) = obj.cnt(I) + 1;
            
            if strcmp(obj.strategy,'UCB1-N') || ...
                    strcmp(obj.strategy.name,'UCB1-V')
                % update sums of sqr rewards
                obj.strategy.q(I) = obj.strategy.q(I) + R^2;
            end
                                      
        end
        
    end
    
    
    
end