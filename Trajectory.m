% Class for holding inputs and states on trajectories.

classdef Trajectory < handle
    
    properties
        
        % time profile
        t
        % number of discretizations
        N
        % discrete set of states on (cts) trajectory
        s
        % nominal u values calculated during trajectory generation
        unom
        % particular algorithm's performance as array of structures
        PERF
        % number of performances attached to a trajectory
        num_per
        % spline points
        sp
        % stopping time
        stoptime
    end
    
    methods (Access = public)
         
        function obj = Trajectory(t,sp,s,unom)            
            obj.t = t;
            obj.N = length(t);
            obj.sp = sp;
            obj.s = s;
            obj.stoptime = length(t);
            obj.unom = unom;
            obj.PERF = struct;
            obj.num_per = 0;
        end
        
        %TODO: simplify this!
        function addPerformance(obj,u,x,x_pred,COST,controller)
            
            if ischar(controller)
                name = controller;
            else
                name = controller.name;
            end
            
            obj.num_per = obj.num_per + 1;
            i = obj.num_per;
            obj.PERF(i).name = name;
            obj.PERF(i).stoptime = size(u,2);
            obj.PERF(i).u = u;
            obj.PERF(i).x = x;
            run = obj.PERF(i).stoptime;
            % add contexts and costs
            % as context, we keep only present state + future state
            x_trim = x(:,1:end-1);
            x_next = x(:,2:end);
            s_next = obj.s(:,2:1+run);
            obj.PERF(i).context = [x_trim;s_next];
            diff = x_next - s_next;
            obj.PERF(i).cost = diag(diff' * COST.Q * diff);
            % nominal model prediction cost
            if ~isempty(x_pred)
                x_pred_next = x_pred(:,2:end);
                diff = x_pred_next - s_next;
                obj.PERF(i).nom_cost = diag(diff' * COST.Q * diff);
            end
            % display SSE error
            sse = sum(obj.PERF(i).cost);
            fprintf('Q-SSE for %s is %f \n', name, sse);
            % add to controller
            if ~ischar(controller)
                controller.record(sse);
            end
        end
                
        % plot the performance of a complete run
        function plot_output(obj,model)
            d = model.dim_tr;
            for i = 1:obj.num_per
                figure;
                title('Desired trj vs. Real trj');
                perf = obj.PERF(i);
                x = perf.x;
                for j = 1:length(d)
                    subplot(2,2,j); 
                    plot(obj.t, obj.s(d(j),:), '-r', obj.t, x(d(j),:), '-g');
                    legend('nominal (desired) trj', 'trj followed');
                end
                subplot(2,2,[3 4]);
                plot(obj.s(d(1),:), obj.s(d(2),:), '-r', x(d(1),:), x(d(2),:), '-g');
                legend('nominal (desired) trj', 'trj followed');
            end
        end
        
        % plot the performance of particular runs
        % till stopping time
        function plot_learning_output(obj,model)
            d = model.dim_tr;
            for i = 2:obj.num_per
                name = obj.PERF(i).name;
                run = obj.PERF(i).stoptime + 1;
                str = sprintf('%s Episode %d: Desired trj vs. Real trj',...
                               name,i-1);
                figure('name', str,'NumberTitle','off');
                perf = obj.PERF(i);
                x = perf.x;
                for j = 1:length(d)
                    subplot(2,2,j); 
                    plot(obj.t, obj.s(d(j),:), '-r', ...
                         obj.t(1:run), x(d(j),:), '-g');
                    legend('nominal (desired) trj', 'trj followed');
                end
                subplot(2,2,[3 4]);
                plot(obj.s(d(1),:), obj.s(d(2),:), '-r', ...
                     x(d(1),:), x(d(2),:), '-g');
                legend('nominal (desired) trj', 'trj followed');
            end
        end
        
        %plot the inputs of a complete run
        function plot_input(obj,model)
            
            dim_u = model.dim_u;
            for i = 1:obj.num_per
                figure;
                perf = obj.PERF(i);
                for j = 1:dim_u;
                    subplot(2,1,j);
                    plot(obj.t(1:end-1),perf.u(j,:),'-r'); 
                    grid on;
                end
                title('Inputs along the trajectory');
            end
        end
        
        %plot the inputs of a particular run
        %till stopping time
        function plot_learning_input(obj,model)
            
            dim_u = model.dim_u;
            for i = 2:obj.num_per
                run = obj.PERF(i).stoptime;
                str = sprintf('Episode %d: Inputs along the trajectory',...
                              i-1);
                figure('name', str,'NumberTitle','off');
                perf = obj.PERF(i);
                for j = 1:dim_u;
                    subplot(2,1,j);
                    plot(obj.t(1:run),perf.u(j,:),'-r'); 
                    grid on;
                end
            end
        end
        
    end
end