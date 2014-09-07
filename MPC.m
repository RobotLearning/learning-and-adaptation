% (Nonlinear) MPC implementing the Controller superclass.
% TODO: very inefficient optimization 

classdef MPC < Controller
    
    
    properties
        
        % number of total episodes so far
        epi
        % color of particular controller
        color
        % name of the particular controller
        name
        % costs incurred (Q-SSE)
        sse
        
        % horizon to consider for the optimization
        horizon
        % dimension of the states
        dim_x
        % dimension of the inputs
        dim_u
        % last input applied
        u_last
        % current t
        tnow
        
    end
    
    methods
    
    function obj = MPC(model,hor)
    
            obj.epi = 0;
            obj.color = 'r';
            obj.name = 'MPC';
            obj.sse = 0;
            
            obj.dim_x = model.dim_x;
            obj.dim_u = model.dim_u;            
            obj.horizon = hor;
            
    end
    
    function obj = update_horizon(obj,trj)
        
            t_now = obj.tnow;
            % do not consider the full horizon towards the end of
            % trajectory
            if t_now > trj.N - obj.horizon
                obj.horizon = trj.N - t_now;
            end
    end
    
    function umin = control(obj,i,trj,model,x)
    
            % fmincon has the advantage of starting from a nonzero value
            if ~isempty(obj.u_last)
                u0 = obj.u_last;
            else
                u0 = trj.unom(:,1);
                %u0 = zeros(model.dim_u,1);
            end
            
            obj.tnow = i;
            u0 = repmat(u0,obj.horizon,1);
            bnd = model.bound;
            bounds = repmat(bnd,obj.horizon,1);
            % update horizon
            obj.update_horizon(trj);
            
            % CALL DIRECT
            %{
            % Send options to Direct 
            options.showits   = 0;
            options.tol       = 0.05;
            options.maxevals  = 200;
            options.maxits    = 100;
            options.maxdeep   = 100;
            % Pass function as part of a Matlab Structure
            problem.f = @(u) obj.rhc(u,trj,model,x);
            [~,umin,hist] = Direct(problem,bounds,options); %#ok
            %disp('DIRECT found xmin at:'); umin
            %figure(2); plot(hist(:,2),hist(:,3),'-*');
            %}

            % CALL FMINCON
            %%{
            f = @(u) obj.rhc(u,trj,model,x);
            %opts = optimset('Display', 'off', 'Algorithm','interior-point');
            opts = optimset('Display', 'off', 'Algorithm','sqp');
            umin = fmincon(f,u0,[],[],[],[],bounds(:,1),bounds(:,2),[],opts);
            % only take the first input at time t
            umin = umin(1:obj.dim_u);
            %disp('FMINCON found xmin at:'); xmin
            %}

            % finish up
            obj.u_last = umin;
        
    end
    
    function val = rhc(obj,u,trj,model,x)

            % get necessary variables out
            dim = obj.dim_u;
            hor = obj.horizon;
            t_now = obj.tnow;

            % minimize the combined cost over the horizon
            nom_err = 0;
            
            x_pred = x;
            for i = 1:hor              
                % use nominal model prediction error
                x_pred = model.predict(t_now,x_pred,...
                                       u((dim*(i-1)+1):(dim*i)));            
                diff = x_pred - trj.s(:,t_now+i);
                nom_err = nom_err + diff' * model.COST.Q * diff;
            end
            val = nom_err;
            % penalize input
            %a = 1.5 * 1e-6;
            %R = a * diag(length(u));
            %pen = u'*R*u;
            %val = nom_err + pen;

    end
    
    % cheats using the "real" dynamics
    function val = cheat(obj,u,trj,model,x)

            % get necessary variables out
            dim = obj.dim_u;
            hor = obj.horizon;
            t_now = obj.tnow;

            % minimize the combined cost over the horizon
            nom_err = 0;
            
            x_pred = x;
            for i = 1:hor              
                % use nominal model prediction error
                x_pred = model.evolve(t_now,x_pred,...
                                      u((dim*(i-1)+1):(dim*i)));            
                diff = x_pred - trj.s(:,t_now+i);
                nom_err = nom_err + diff' * model.COST.Q * diff;
            end
            val = nom_err;
            % penalize input
            %a = 1.5 * 1e-6;
            %R = a * diag(length(u));
            %pen = u'*R*u;
            %val = nom_err + pen;

    end
        
    end
end