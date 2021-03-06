% Model superclass acts as an interface
% for holding specific submodel classes.

classdef (Abstract) Model < handle
    
    properties (Abstract)
        
        % parameters structure
        PAR
        % constraints structure
        CON
        % disturbance structure
        DIST
        % rectangular input bound (matrix)
        bound
        % cost function structure (handle and weight matrix Q)
        COST
        % flag for previously saved experiments
        flg_exp
        % file for saving trajectories 
        expSavefile
        % flag for previously saved estimation results
        flg_est
        % file for saving estimation results
        estSavefile
    end
    
    properties
        % simulation time step
        h
        % dimension of the x vector
        dim_x
        % dimension of the action space
        dim_u
        % 2 dimensions, y and z, are tracked only
        dim_tr
        % noise variance of the model (white noise)
        eps
    end
    
    % methods to be implemented
    methods (Abstract)
        
        % generate trajectories
        trajectory(t,sp,s,unom)
        % plot trajectory tracking results
        plot(obj,trjs)
        % set a nominal model
        nominal(t,obj,x,u)
        % set a disturbance model
        disturbance(t,obj,x,u)
        % linearize around trajectory
        linearize(obj,trj)
        % get lifted constraints
        lift(obj,trj)
        
    end
    
    % methods that can be implemented here in abstract class
    methods (Access = public)
        
        % simulates whole trajectory
        function x = simulate(obj,t,x0,us,fun)
            N = length(t)-1;
            h = t(2)-t(1);
            x = zeros(length(x0),N+1);
            x(:,1) = x0;
            for i = 1:N
                x(:,i+1) = step(obj,t(i),x(:,i),us(:,i),fun);
                % no constraint checking
            end
        end
        
        % one step simulation along the trajectory
        % TODO: is it correct to assume u constant?
        function next = step(obj,t,prev,u,fun)
            
            % get trajectory of states
            % using classical Runge-Kutta method (RK4)
            k1 = obj.h * fun(t,obj,prev,u);
            x_k1 = prev + k1/2;
            k2 = obj.h * fun(t,obj,x_k1,u);
            x_k2 = prev + k2/2;
            k3 = obj.h * fun(t,obj,x_k2,u);
            x_k3 = prev + k3;
            k4 = obj.h * fun(t,obj,x_k3,u);
            next = prev + (k1 + 2*k2 + 2*k3 + k4)/6;

        end
        
        % predict next states using function
        % useful to calculate nominal model prediction error
        function x_pre = predict_full(obj,t,x,us,fun)
            N = length(t)-1;
            x_pre = zeros(size(x,1),N+1);
            x_pre(:,1) = x(:,1);
            for i = 1:N
                x_pre(:,i+1) = step(obj,t(i),x(:,i),us(:,i),fun);
            end            
        end
        
        % useful to propagate one full iteration of
        % ILC input sequence
        function x_next = evolve_full(obj,t,x,us)
            fun = @(t,obj,x,u) nominal(t,obj,x,u,false) ...
                             + disturbance(t,obj,x,u);
            x0 = x(:,1);
            x_next = simulate(obj,t,x0,us,fun);
        end
        
        % predict one step using nominal model
        function x_pre = predict(obj,t,x,u)
            fun = @(t,obj,x,u) nominal(t,obj,x,u,false);            
            x_pre = step(obj,t,x,u,fun);
        end
        
        % evolve one step
        function x_next = evolve(obj,t,x,u)            
            fun = @(t,obj,x,u) nominal(t,obj,x,u,false) ...
                             + disturbance(t,obj,x,u);
            x_next = step(obj,t,x,u,fun);
        end
        
        % generate examples for hyperparameter estimation
        % and then estimate with maximum likelihood
        function est = experiment(obj,shapes,coords)
           if obj.flg_exp 
               % apply function to each cell in cell array
               trjs = cellfun(@(shape,coord) trajectory(obj,shape,coord), ...
                                shapes, coords, 'UniformOutput', false);
               % plot trajectories
               obj.plot(trjs);
               % save trajectories cell
               save(obj.expSavefile,'trjs');               
           else
               % load trajectories
               load(obj.expSavefile);
           end
           
           if obj.flg_est
               est = Estimator(obj);
               [x,y] = est.processCtrlExp(trjs);
               est.mle(x',y);
               % save estimation results
               save(obj.estSavefile,'est');
           else
               % load hyperparameters
               load(obj.estSavefile);
           end
        end       
        
        
    end
end