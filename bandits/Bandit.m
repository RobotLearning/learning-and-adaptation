%% Bandit class for implementing different strategies

classdef Bandit < handle
    
    properties
       
        % strategies for exploring/exploiting
        strategy
        % mean and variance statistics for each arm
        mean
        % counters for number of plays so far
        cnt
        
    end
    
    methods
        
        %% Initialize Bandit with a strategy
        function obj = Bandit(arms,str)
           
            obj.mean = zeros(arms,1);
            obj.strategy = str;
            obj.cnt = zeros(arms,1);
            
            if strcmp(obj.strategy.name,'UCBV')
                obj.strategy.q = zeros(arms,1); % sum of sqr rewards
            end
            
            if strcmp(obj.strategy.name,'UCB')
                % fix the range of rewards
                obj.strategy.b = 1;
            end
            
            if strcmp(obj.strategy.name,'UCBV-Regularized')
               obj.strategy.q = zeros(arms,1); % sum of sqr rewards
               obj.strategy.lambda = str.lambda;
               obj.strategy.last_arm = 1;
            end
            
            % prepare the priors
            if strfind(obj.strategy.name,'Thompson')
                obj.cnt = ones(arms,1);
                obj.strategy.prior.var.a = obj.strategy.a * ones(arms,1);
                obj.strategy.prior.var.b = obj.strategy.b * ones(arms,1);
            end
            
            if strfind(obj.strategy.name,'Cautious')
                % add sampling number
                obj.strategy.sample_num = obj.strategy.sample_num;
                obj.strategy.last_arm = 1;
            end
            
            if strcmp(obj.strategy.name,'Thompson-Regularized')
                obj.strategy.lambda = str.lambda;
                obj.strategy.last_arm = 1;
            end
            
        end
        
        %% Strategies for playing
        function idx = play(obj)
            
            n_arms = length(obj.mean);
            switch obj.strategy.name
                case 'EPS-GREEDY'
                    idx = obj.epsGreedy(n_arms);
                case 'EPS-FIRST'
                    idx = obj.epsFirst(n_arms);
                case 'UCB'
                    idx = obj.ucb1(n_arms);
                case 'UCBV'
                    idx = obj.ucb1Variance(n_arms);
                case 'Thompson-Normal'
                    idx = obj.thompson_gauss(n_arms);
                case 'Thompson-Cautious'
                    idx = obj.thompson_cautious(n_arms);
                case 'Thompson-Regularized'
                    idx = obj.thompson_regular(n_arms);
                case 'UCBV-Regularized'
                    idx = obj.ucb_vregular(n_arms);
                otherwise
                    error('Alg not implemented!');
            end
            
        end
        
        % Thompson sampling for Gaussian distributions
        function idx = thompson_gauss(obj,arms)
            % just sample from the posterior which 
            % has been updated from the prior
            
            % first sample from the variance distributions
            a = obj.strategy.prior.var.a;
            b = obj.strategy.prior.var.b;
            vars = zeros(arms,1);
            for i = 1:arms
                vars(i) = inv_gamma(a(i),b(i)); % TODO:
            end
            means = obj.mean + sqrt(vars./obj.cnt).*randn(arms,1);
             
             [~,idx] = max(means);
            
        end
        
        % Thompson sampling regularized w.r.t. last-arm
        function idx = thompson_regular(obj,n_arms)

            lambda = obj.strategy.lambda;
            last_arm = obj.strategy.last_arm;
            % first sample from the variance distributions
            a = obj.strategy.prior.var.a;
            b = obj.strategy.prior.var.b;
            vars = zeros(n_arms,1);
            for i = 1:n_arms
                vars(i) = inv_gamma(a(i),b(i)); % TODO:
            end
            arms = 1:n_arms;
            cost = lambda * abs(arms - last_arm);
            means = obj.mean - cost(:) + sqrt(vars./obj.cnt).*randn(n_arms,1);
            [~,idx] = max(means);
            obj.strategy.last_arm = idx;
            
        end
        
        % Thompson sampling with resampling
        function idx_closest = thompson_cautious(obj,arms)
            % just sample from the posterior which 
            % has been updated from the prior
            
            t = sum(obj.cnt) + 1;
            % first sample from the variance distributions
            N = obj.strategy.sample_num(t);
            last_arm = obj.strategy.last_arm;
            a = obj.strategy.prior.var.a;
            b = obj.strategy.prior.var.b;
            vars = zeros(arms,1);
            idx = zeros(N,1);            
            
            for n = 1:N % repeat sampling from maximum
                for i = 1:arms
                    vars(i) = inv_gamma(a(i),b(i)); 
                end
                means = obj.mean + sqrt(vars./obj.cnt).*randn(arms,1);
                [~,idx(n)] = max(means);
            end
            % choose closest index to last_arm
            [~,I] = min(abs(idx - last_arm)); % find closest
            % last arm
            obj.strategy.last_arm = idx(I);
            idx_closest = idx(I);
        end
        
        % regularized ucb strategy
        function idx = ucb_vregular(obj,n_arms)
            
            arms = 1:n_arms;
            lambda = obj.strategy.lambda;
            last_arm = obj.strategy.last_arm;
            t = sum(obj.cnt) + 1;
            % for variance aware ucb
            qs = obj.strategy.q ./ obj.cnt;
            n = obj.cnt - 1;
            n(n == 0) = 1;
            V = (qs - obj.mean.^2)./n;
            rho = 1;
            if t <= n_arms
                idx = t;
            else
                cost = lambda * abs(arms - last_arm);
                biasEst = sqrt(rho*log(t).*V)./sqrt(obj.cnt);
                [~,idx] = max(obj.mean - cost(:) + biasEst);
            end
            obj.strategy.last_arm = idx;
                
        end
        
        % classical ucb strategy
        function idx = ucb1(obj,arms)
            
            % estimate the range of the rewards
            b = obj.strategy.b;
            t = sum(obj.cnt) + 1;
            % make sure each machine is played once
            rho = obj.strategy.rho;
            if t <= arms
                idx = t;
            else
                [~,idx] = max(obj.mean + sqrt(b^2*rho.*log(t))./sqrt(obj.cnt));
            end
                
        end
        
        % classical ucb with variance estimate
        function idx = ucb1Variance(obj,arms)
            
            t = sum(obj.cnt) + 1;
            % for variance aware ucb
            qs = obj.strategy.q ./ obj.cnt;
            n = obj.cnt - 1;
            n(n == 0) = 1;
            V = (qs - obj.mean.^2)./n;
            rho = 1;
            if t <= arms
                idx = t;
            else
                biasEst = sqrt(rho*log(t).*V)./sqrt(obj.cnt);
                [~,idx] = max(obj.mean + biasEst);
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
            
            if strfind(obj.strategy.name,'Thompson')
                % get the posteriors
                obj.strategy.prior.var.a(I) = obj.strategy.prior.var.a(I) + ...
                                           1/2;
                obj.strategy.prior.var.b(I) = obj.strategy.prior.var.b(I) + ...
               0.5 * (obj.cnt(I) / (obj.cnt(I) + 1)) * (obj.mean(I) - R)^2;
            end
            
            % update mean and add one to chosen arm's counter
            obj.mean(I) = (obj.mean(I)*obj.cnt(I) + R)/(obj.cnt(I)+1);
            obj.cnt(I) = obj.cnt(I) + 1;
            
            if strcmp(obj.strategy.name,'UCBV')
                % update sums of sqr rewards
                obj.strategy.q(I) = obj.strategy.q(I) + R^2;
            end
            
            % update if range seems to be bigger
            if strcmp(obj.strategy,'UCB')
                rangeEst = R;
                obj.strategy.b = max(1,rangeEst);
            end
            
        end
        
    end
    
    
    
end