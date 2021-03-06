%% Bandit class for implementing different strategies

classdef BO < handle
    
    properties
       
        % mesh of points to predict the function at
        mesh
        % strategies for exploring/exploiting
        strategy
        % gaussian process
        gp
        % counters for number of plays so far
        t
        
    end
    
    methods
        
        %% Initialize Bayesian Optimization with a strategy
        function obj = BO(mesh,hp,str)
           
            obj.mesh = mesh;
            obj.gp = GP(hp,[],[]);
            obj.strategy = str;
            obj.t = 1;
            
            if strcmp(obj.strategy.name,'GP-UCB')
                % hoeffding prob. constant
                obj.strategy.delta = str.delta;
            elseif strcmp(obj.strategy.name,'GP-MI')
                obj.strategy.alpha = str.alpha;
                obj.strategy.gamma = 0.0;
            elseif strcmp(obj.strategy.name,'Thompson-Cautious')
                obj.strategy.buffer = str.buffer;
                obj.strategy.last_idx = 1;
            elseif strcmp(obj.strategy.name,'Thompson-Regular')
                obj.strategy.lambda = str.lambda;
                obj.strategy.last_idx = 1;
            elseif strcmp(obj.strategy.name,'Thompson-Plan')
                obj.strategy.buffer = str.buffer;
                obj.strategy.last_arms = ones(1,str.buffer);
                obj.strategy.current_idx = str.buffer;
                % number of times to sample
                %error('Alg not implemented!');
            end
        end
        
        %% Strategies for playing
        function [x,idx] = play(obj)
            
            switch obj.strategy.name
                case 'GP-UCB'
                    [x,idx] = obj.gp_ucb();
                case 'GP-MI'
                    [x,idx] = obj.gp_mi();
                case 'EI'
                    [x,idx] = obj.ei();
                case 'Thompson-Normal'
                    [x,idx] = obj.thompson();
                case 'Thompson-Cautious'
                    [x,idx] = obj.thompson_cautious();
                case 'Thompson-Regular'
                    [x,idx] = obj.thompson_regular();
                case 'Thompson-Plan'
                    [x,idx] = obj.thompson_plan();
                otherwise
                    error('Alg not implemented!');
            end            
        end
        
        function [x,idx] = gp_ucb(obj)
        
           meshsize = length(obj.mesh);
           delta = obj.strategy.delta;
           beta_thm1 = @(t,D,delta) 0.1*log((t^2)*D*(pi^2)/(6*delta));
           beta = beta_thm1(obj.t,meshsize,delta);
           [mu,Sigma] = obj.gp.predict_mesh(obj.mesh);
           s2 = diag(Sigma);
           acq = @(m,s2) m + sqrt(beta * s2);
           [~,idx] = max(acq(mu(:),s2(:)));
           idx = idx(1);
           x = obj.mesh(:,idx);
        end
        
        function [x,idx] = gp_mi(obj)
            alpha = obj.strategy.alpha;
            gamma = obj.strategy.gamma;
            acq = @(m,s2) m + sqrt(alpha * (s2 + gamma)) ...
                            - sqrt(alpha * gamma);
            [mu,Sigma] = obj.gp.predict_mesh(obj.mesh);
            s2 = diag(Sigma);
            [~,idx] = max(acq(mu(:),s2(:)));
            idx = idx(1);
            x = obj.mesh(:,idx);
            obj.strategy.gamma = gamma + s2(idx);
        end
        
        function [x,idx] = ei(obj)
           
            density = @(x) (1/sqrt(2*pi)) .* exp((-1/2).*(x.^2));
            distr = @(x) (1 - erf(x./sqrt(2)))./2;
            acq = @(m,s2) (m - max(m)) .* distr((max(m) - m)./sqrt(s2)) + ...
                      sqrt(s2) .* density((max(m) - m)./sqrt(s2));
            [mu,Sigma] = obj.gp.predict_mesh(obj.mesh);
            s2 = diag(Sigma);
            [~,idx] = max(acq(mu(:),s2(:)));
            idx = idx(1);
            x = obj.mesh(:,idx);            
        end
        
        %% THOMPSON SAMPLING METHODS
        function [x,idx] = thompson(obj)
            meshsize = length(obj.mesh);
            [mu,Sigma] = obj.gp.predict_mesh(obj.mesh);
            % better for too smooth kernels
            [U,S] = eig(Sigma);
            f = mu(:) + 0.2 * U * sqrt(max(0,real(S))) * randn(meshsize,1);
            [~,idx] = max(f);
            idx = idx(1);
            x = obj.mesh(:,idx);
        end
        
        % regularized thompson
        function [x,idx] = thompson_regular(obj)
            last_idx = obj.strategy.last_idx;
            lambda = obj.strategy.lambda;
            meshsize = length(obj.mesh);
            idxs = 1:meshsize;
            [mu,Sigma] = obj.gp.predict_mesh(obj.mesh);
            % better for too smooth kernels
            [U,S] = eig(Sigma);
            dist = obj.mesh(:,idxs) - obj.mesh(:,last_idx);
            cost = sqrt(sum(dist.*dist,1));
            f = mu(:) - lambda * cost(:) + 0.2 * U * sqrt(max(0,real(S))) * randn(meshsize,1);
            [~,idx] = max(f);
            idx = idx(1);
            obj.strategy.last_idx = idx;
            x = obj.mesh(:,idx);
        end        
        
        % cautious thompson
        function [x,idx] = thompson_cautious(obj)
            meshsize = length(obj.mesh);
            buffer = obj.strategy.buffer(obj.t);
            last_idx = obj.strategy.last_idx;
            [mu,Sigma] = obj.gp.predict_mesh(obj.mesh);
            % better for too smooth kernels
            [U,S] = eig(Sigma);
            M = 0.2 * U * sqrt(max(0,real(S)));
            f = zeros(meshsize,buffer);
            for i = 1:buffer
                f(:,i) = mu(:) + M * randn(meshsize,1);
            end
            [~,idx] = max(f);
            % choose closest index to last index taken
            dists = obj.mesh(:,idx) - obj.mesh(:,last_idx);
            [~,I] = min(sum(dists.*dists,1)); 
            % find closest
            idx = idx(I);
            obj.strategy.last_idx = idx;
            x = obj.mesh(:,idx);
        end        
        
        % planning thompson
         function [x,idx] = thompson_plan(obj)
            meshsize = length(obj.mesh);
            buffer = obj.strategy.buffer;
            last_arms = obj.strategy.last_arms;
            current_idx = obj.strategy.current_idx;
            if current_idx == buffer
                obj.strategy.current_idx = 1;
                [mu,Sigma] = obj.gp.predict_mesh(obj.mesh);
                % better for too smooth kernels
                [U,S] = eig(Sigma);
                M = U * sqrt(max(0,real(S)));
                f = zeros(meshsize,buffer);
                for i = 1:buffer
                    f(:,i) = mu(:) + M * randn(meshsize,1);
                end
                [~,idx] = max(f);
                % FIND NEAREST NEIGHBORS
                % choose closest index to last index taken
                list = idx;
                last_idx = current_idx;
                for i = 1:buffer
                    dists = obj.mesh(:,list) - obj.mesh(:,last_idx);
                    % remove an element 
                    [~,I] = min(sum(dists.*dists,1));
                    obj.strategy.last_arms(i) = list(I);
                    last_idx = list(I);
                    list = setdiff(list,list(I));
                end
            else
                obj.strategy.current_idx = current_idx + 1;
                x = obj.mesh(:,last_arms(current_idx+1));
                idx = last_arms(current_idx+1);
            end
         end
        
        
        %% Get rewards and update GP
        function update(obj,x,y)
            
            % update GP and other stats here
            % update mean and variance of GP
            obj.t = obj.t + 1;
            obj.gp.update(x,y);
            
        end
        
    end
    
    
    
end