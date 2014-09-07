% Trajectory Tracking using Gaussian Process Optimization

classdef TGP < Controller
    
    properties (Constant)
        
        % construct beta_t using Thm 2 in A.Krause's paper
        % Hoeffding inequality constant
        delta = 0.1;
        % these constants encode smoothness of functions
        % drawn from the kernel (see GPUCB paper)
        a = 1; b = 1;
        % Beta multiplier
        % more aggression necessary for tackling real problems
        beta_mult = 0.1;
        % Stopping time
        epsilon = 0.1;
        
    end
    
    properties
        
        % flag for using Rasmussen's gp inference or custom gp inference
        flg_gp
        % in custom mode use GP structure instead of estimator class
        GP
        % number of total episodes so far
        epi
        % color of particular controller
        color
        % name of the particular controller
        name
        % costs incurred (Q-SSE)
        sse
        
        % holds the past cost differences
        ys
        % holds the past contexts
        ctxs
        % holds the past control inputs
        us
        % estimator class
        est
        % next state to go
        s_next
        % current position
        x_now
        % current time
        t_now
        % total time stage throughout all past+current episodes
        T
    end
    
    properties
        
        % Bounds for compact context-action space
        r
        % Dimension of context-action space
        d
        % Handle for scaling exploration-exploitation
        beta_f = @(mult,t,delta,d,a,b,r) mult * ...
                  (2*log((t^2)*2*(pi^2)/(3*delta)) +...
                   2*d*log((t^2)*d*b*r*sqrt(log(4*d*a/delta))));        
    end
    
    methods
        
        function obj = TGP(model,est,flg_inf)
            
            obj.T = 1;
            obj.epi = 0;
            obj.color = 'b';
            obj.name = 'TGP';
            obj.sse = 0;
            obj.flg_gp = flg_inf;
            
            dim_ctx = 2 * model.dim_x;
            dim_u = model.dim_u;
            
            obj.r = max(max(model.bound));
            obj.d = dim_ctx + dim_u;
            obj.est = est;
            obj.ys = double.empty(0,1);
            obj.ctxs = double.empty(dim_ctx,0);
            obj.us = double.empty(dim_u,0);
            
            if flg_inf
                obj.GP = struct();
                obj.custom(model);
            end
        end
        
        % prepare to use custom gp-inference
        % TODO: automate this properly, still using Rasmussen's
        % hp-estimation procedure
        function custom(obj,model)
            
            % release hyperparameters
            dim_ctx = 2 * model.dim_x;
            hp_ctx = Inf(1,dim_ctx);
            maskCtx = obj.est.cov{2}{1}{2}{1};
            ctx_l = sum(maskCtx);
            hp_ctx(maskCtx > 0) = obj.est.hp.cov(1:ctx_l);
            hp_act = [obj.est.hp.cov(end),Inf];

            % variance of the estimated noise
            var_n = exp(obj.est.hp.lik)^2;
            % kernel function for the contexts
            ctx_hp = exp(hp_ctx);
            ctx_scale = exp(obj.est.hp.cov(ctx_l+1))^2;
            ctx_type = 'squared exponential ard';
            ctx_kern = @(x1,x2) ctx_scale * kernel(x1(1:dim_ctx), ...
                                       x2(1:dim_ctx), ctx_hp, ctx_type);
            % kernel function for the actions
            act_hp = exp(hp_act);
            act_scale = 1;
            act_type = 'linear ard';
            act_kern = @(x1,x2) act_scale * kernel(x1(dim_ctx+1:end), ...
                                       x2(dim_ctx+1:end), ...
                                       act_hp, act_type);
            kern = @(x1,x2) act_kern(x1,x2) * ctx_kern(x1,x2);
            % previous inverse
            %Kinv = 1/(1+exp(est.hp.lik)^2);
            K = 0;
            
            %obj.GP.Kinv = Kinv;
            obj.GP.K = K;
            obj.GP.kern = kern;
            obj.GP.hp.act = hp_act;
            obj.GP.hp.ctx = hp_ctx;
            obj.GP.hp.var_n = var_n;
        end
        
        % Scale exploration-exploitation at a particular time t
        function val = beta(obj,i)
            val = obj.beta_f(obj.beta_mult,i,obj.delta,obj.d,obj.a,...
                             obj.b,obj.r);%
            
        end
        
        % update past contexts
        function update_ctx(obj)
            
            obj.T = obj.T + 1;
            ctx = [obj.x_now;obj.s_next];
            obj.ctxs = [obj.ctxs, ctx];            
        end
        
        % update past inputs
        function update_input(obj,u)
            obj.us = [obj.us, u];
        end
        
        % update past costs
        function update_data(obj,model,cost,u)
            
            % subtract nominal model prediction error
            x_pred = model.predict(obj.t_now,obj.x_now,u);  
            diff = x_pred - obj.s_next;
            nom_cost = diff' * model.COST.Q * diff;
            y = cost - nom_cost;
            obj.ys = [obj.ys; y];
        end
            
        
        function u = control(obj,i,trj,model,x)
            obj.s_next = trj.s(:,i+1);
            obj.t_now = trj.t(:,i);
            obj.x_now = x;
            beta = obj.beta(obj.T);
            
            % fmincon has the advantage of starting from a nonzero value
            if ~isempty(obj.us)
                u0 = obj.us(:,end); 
            else
                u0 = trj.unom(:,1);
                %u0 = zeros(model.dim_u,1);
            end
            
            % custom mode
            if obj.flg_gp
               K = obj.GP.K;
               kern = obj.GP.kern;
               var_noise = obj.GP.hp.var_n;
               cxu = [obj.ctxs; obj.us];
               if obj.T > 1
                   mat = ker_matrix_iter(obj.T-1,cxu,kern,K);
                   % add noise to last point
                   mat(end,end) = mat(end,end) + var_noise;
                   obj.GP.K = mat;
               end
               f = @(x) obj.gpucb(beta,model,x);
            else
               f = @(x) obj.quantile(beta,model,x);
            end
            
            bnds = model.bound;
            % CALL FMINCON
            %%{

            opts = optimset('Display', 'off', 'Algorithm','interior-point');
            %opts = optimset('Display', 'off', 'Algorithm','sqp');
            %disp('FMINCON found xmin at:'); xmin
            u = fmincon(f,u0,[],[],[],[],bnds(:,1),bnds(:,2),[],opts);
            %}
            
            % finish up
            obj.update_ctx();
            obj.update_input(u);
        end
        
        % calculates mean - sqrt(beta) * sigma
        % pt is the test point
        %
        % this is called by the control function many times 
        % so it should not spend idle time processing info
        function val = quantile(obj,beta,model,u)
            
            % add nominal model prediction error to gp
            x_pred = model.predict(obj.t_now,obj.x_now,u);            
            diff = x_pred - obj.s_next;
            nom_cost = diff' * model.COST.Q * diff;
            
            % call gp
            ctx = [obj.x_now;obj.s_next];
            pt = [ctx; u]';
            xs = [obj.ctxs; obj.us]';
            e = obj.est;
            [m, s2] = gp(e.hp, e.inf, e.mean, e.cov, e.lik, xs, obj.ys, pt);
            val = nom_cost + m - sqrt(beta * s2);
        end
        
        % custom gp inference
        function val = gpucb(obj,beta,model,u)
            
            % add nominal model prediction error to gp
            x_pred = model.predict(obj.t_now,obj.x_now,u);            
            diff = x_pred - obj.s_next;
            nom_cost = diff' * model.COST.Q * diff;

            % call gp after getting one cost (for conditioning)
            if obj.T > 1
                ucb_val = obj.mu_beta_sigma(beta,u);
                val = nom_cost + ucb_val;
            else
                val = nom_cost;
            end

            % debug beta level (exploration)
            %fprintf('nom = %f, m = %f, s2 = %f, explore = %f.\n', ...
            %        nom_cost, m, s2, sqrt(beta * s2));

        end

        % custom ucb function
        function val = mu_beta_sigma(obj,beta,u)

            % get necessary variables out
            i = obj.T - 1;
            ctx = [obj.x_now;obj.s_next];
            kern = obj.GP.kern;
            %mat = obj.GP.Kinv;
            mat = obj.GP.K;
            cxu = [obj.ctxs; obj.us];
            x = [ctx; u]';
            costs = obj.ys;
            
            % calculate kernel vector k and store
            vec = ker_vector(i,x,cxu,kern);
            % create mu function 
            mu = vec' * (mat \ costs(:));
            %mu = vec' * mat * costs(:);
            % create sigma function
            s2 = max(kern(x,x) - vec' * (mat \ vec), 0);
            %s2 = max(kern(x,x) - vec' * mat * vec, 0);
            % pass mu and sigma functions as function handles
            val = mu - sqrt(beta * s2);

        end
        
    end
    
end