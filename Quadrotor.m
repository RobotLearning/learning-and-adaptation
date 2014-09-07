% Quadrotor model implementing the Model superclass

classdef Quadrotor < Model

    properties   
        % parameters structure
        PAR
        % constraints structure
        CON
        % input bound
        bound
        % disturbance structure
        DIST
        % cost function structure (handle and weight matrix)
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
    
    methods
        
        % copies the parameter values inside the structure
        function set.PAR(obj, STR)  
            % initalize all fields to 0
            obj.PAR = struct('Iy',0,'L',0,'m',0,'g',0);
            obj.PAR.Quad = struct('A',0);
                         
            % check that the input has all the fields
            % TODO: is there a better way?
            assert(all(strcmp(fieldnames(obj.PAR), fieldnames(STR))));
            obj.PAR = STR;
        end
        
        % copies the constraint values inside the structure
        function set.CON(obj, STR)
            % initialize all fields to 0
            obj.CON = struct('fi_max',0,'fi_min',0,'fi_dot_max',0,...
                             'fmax',0,'fmin',0,'f_dot_max',0, 'phi_max',0,...
                             'phi_dot_max',0,'phi_ddot_max',0,...
                             'NumT',0);
            % check that the input has all the fields
            assert(all(strcmp(fieldnames(obj.CON), fieldnames(obj.CON))));
            obj.CON = STR;
            
        end
        
        % copies the rectangular bounds on R2 input space
        function set.bound(obj, mat)
            
            if isempty(obj.bound)
                obj.bound = Inf * [-1,1;
                                   -1,1];
            end
            
            % make sure the values are tighter
            obj.bound(:,1) = max(obj.bound(:,1), mat(:,1));
            obj.bound(:,2) = min(obj.bound(:,2), mat(:,2));
            
        end
        
        % change flag for saved experiment
        function set.flg_exp(obj,flg)
            obj.flg_exp = flg;
        end
        
        % change flag for saved estimation results
        function set.flg_est(obj,flg)
            obj.flg_est = flg;
        end
        
        % change the cost function
        function set.COST(obj, Q)
            obj.COST.Q = Q;
            obj.COST.fnc = @(x1,x2) (x1-x2)'*Q*(x1-x2);
        end
        
        % change the disturbance type
        function set.DIST(obj, STR)
            allowable_dist = {'gravity','wind','actuator'};
            assert(any(strcmp(STR.Type,allowable_dist)));
            % do not compare fields here
            obj.DIST = STR;
        end
    end
    
    methods
        
        % constructor for convenience
        % TODO: divide into several methods?
        function obj = Quadrotor(par,dist,con,Q)
            
            % simulation time step
            obj.h = 0.02;
            % dimension of the x vector
            obj.dim_x = 5;
            % dimension of the action space
            obj.dim_u = 2;
            % 2 dimensions, y and z, are tracked only
            obj.dim_tr = [1,3]; 
            % noise variance of the model (white noise)
            obj.eps = 0.00003;
            
            % set object parameter
            obj.PAR = par;
            % set object disturbance
            obj.DIST = dist;
            % set object constraints
            obj.CON = con;
            
            % bounds of the inputs
            bnd = [con.fmin, con.fmax; 
                  -con.phi_dot_max, con.phi_dot_max]; 
            % set object bounds
            obj.bound = bnd;
                     
            % cost function handle
            obj.COST = Q;
            
            % experiment file
            obj.expSavefile = 'quad_exp.mat';
            obj.flg_exp = false;
            
            if ~exist(obj.expSavefile, 'file')
                % no experiments have been performed so far
                obj.flg_exp = true;
            end
            
            % estimation file
            obj.estSavefile = 'quad_est.mat';
            obj.flg_est = false;
            
            if ~exist(obj.estSavefile, 'file')
                % no estimation has been performed so far
                obj.flg_est = true;
            end
        end
        
        % provides nominal model
        function [x_dot,varargout] = nominal(~,obj,x,u,flg)
            % x_dot = quadrocopterNominalDynamics(t,x,u)
            % differential equation of the quadrocopter
            % Dynamics
            % x_dot = Ax + B(x)u + C
            A = [0 1 0 0 0;
                 zeros(1,5);
                 0 0 0 1 0;
                 zeros(2,5)];
            B = [0          0;
                 -sin(x(5)) 0;
                 0          0;
                 cos(x(5))  0;
                 0          1];
            C = [0; 0; 0; -obj.PAR.g; 0];
            x_dot = A*x + B*u + C;
            
            % return jacobian matrices
            if flg
                vec = [0; -u(1)*cos(x(5)); 0; u(1)*sin(x(5)); 0];
                dfdx = [A(:,1:4), vec];
                varargout{1} = dfdx;
                varargout{2} = B;
            end
        end
        
        % provides disturbance to be added to the nominal model
        function delta = disturbance(~,obj,x,u)
             
            switch obj.DIST.Type
                case 'gravity'
                    F_z = obj.PAR.g - obj.DIST.g;
                    delta = [0; 0; 0; F_z; 0];                
                case 'wind'
                    F_y = obj.DIST.BlowPressure * obj.PAR.Quad.A * ...
                          sin(obj.DIST.Angle + x(5)) * cos(obj.DIST.Angle);
                    F_z = obj.DIST.BlowPressure * obj.PAR.Quad.A * ...
                          sin(obj.DIST.Angle + x(5)) * sin(obj.DIST.Angle);
                    delta = [0; F_y; 0; F_z; 0];
                case 'actuator'
                    if ~isfield(obj.DIST, 'k')
                        obj.DIST.k = 0.1; % multiply B matrix by 1+k
                    end
                    delta_B = [0 0;
                              -obj.DIST.k*sin(x(5)) 0;
                               0 0;
                               obj.DIST.k*cos(x(5)) 0;
                               0 1];
                    delta = delta_B * u;
            end
        end
        
        % linearizes the nominal dynamics around the trajectory
        function [A,B] = linearize(obj,trj)
            
            N = trj.N - 1; 
            t = trj.t;
            s = trj.s;
            dimx = obj.dim_x;
            dimu = obj.dim_u;
            unom = trj.unom;
            A = zeros(dimx,dimx,N);
            B = zeros(dimx,dimu,N);
            for i = 1:N
                [~,A(:,:,i), B(:,:,i)] = nominal(t(i),obj,s(:,i),...
                                                 unom(:,i),true);
                % get discrete approximation from jacobian
                % crude approximation
                %A(:,:,i) = eye(dimx,dimx) + obj.h * A(:,:,i);
                %B(:,i) = obj.h * B(:,i);
                % exact matrix calculation 
                Mat = [A(:,:,i), B(:,:,i); zeros(dimu, dimx + dimu)];
                MD = expm(obj.h * Mat);
                A(:,:,i) = MD(1:dimx,1:dimx);
                B(:,:,i) = MD(1:dimx,dimx+1:end);
            end
        end
        
        % get lifted model constraints
        function [umin,umax,L,q] = lift(obj,trj)
            
            N = trj.N - 1; 
            u_trj = trj.unom(:,1:N);
            %dimx = obj.dim_x;
            dimu = obj.dim_u;
            umin(1,:) = obj.CON.fmin - u_trj(1,:);
            umin(2,:) = -obj.CON.phi_dot_max - u_trj(2,:);
            umax(1,:) = obj.CON.fmax - u_trj(1,:);
            umax(2,:) = obj.CON.phi_dot_max - u_trj(2,:);

            % arrange them in a format suitable for optimization
            umin = umin(:);
            umax = umax(:);
            
            % construct D
            D = (diag(ones(1,dimu*(N-1)),dimu) - eye(dimu*N))/obj.h;
            D = D(1:end-dimu,:);
            % construct L1 and L2
            L1 = zeros(N-1, N*dimu);
            L2 = zeros(N-1, N*dimu);
            a = 1/4;
            b = obj.PAR.Iy/(2 * obj.PAR.m * obj.PAR.L); 
            b = b/obj.h; %b_bar
            vec1 = [a -b 0 b];
            vec2 = [a b 0 -b];
            for i = 1:N-1
                L1(i,:) = [zeros(1,(i-1)*dimu), vec1, zeros(1,(N-i-1)*dimu)];
                L2(i,:) = [zeros(1,(i-1)*dimu), vec2, zeros(1,(N-i-1)*dimu)];
            end
            u_dot_max = [4*obj.CON.fi_dot_max; obj.CON.phi_ddot_max];
            U_dot_max = repmat(u_dot_max,N-1,1);
            u_star = u_trj(:); 

            L = [D; -D; L1; -L1; L2; -L2];
            q = [U_dot_max - D*u_star; 
                 U_dot_max + D*u_star;
                 obj.CON.fmax - L1*u_star;
                 -obj.CON.fmin + L1*u_star;
                 obj.CON.fmax - L2*u_star;
                 -obj.CON.fmin + L2*u_star];    
            
        end

        % wrapper for the splines-based trajectory generator method
        % using the differential flatness of the dynamics
        function Traj = trajectory(obj,shape,spline)

            spline_x = spline(1,:);
            spline_y = spline(2,:);
            [t,u_nom,x_nom] = quad_traj_gen(shape,obj.PAR,obj.CON,...
                                            obj.h,spline_x,spline_y);
            fun0 = @(t,obj,x,u) nominal(t,obj,x,u,false);
            fun = @(t,obj,x,u) nominal(t,obj,x,u,false) ...
                             + disturbance(t,obj,x,u);
            x_real = simulate(obj,t,x_nom(:,1),u_nom,fun);
            x_pred = predict_full(obj,t,x_real,u_nom,fun0);
            Traj = Trajectory(t,spline,x_nom,u_nom);
            cost = obj.COST;
            % trajectory generator generates an extra input
            Traj.addPerformance(u_nom(:,1:end-1),x_real,...
                                x_pred,cost,'splines');
            % TODO add noise to x! (as observed variable)
        end
        
        % plot trajectories
        function plot(obj,trjs)
            
            n = length(trjs);
            for i = 1:n
                trjs{i}.plot_output(obj);
                trjs{i}.plot_input(obj);
            end                    
        end
        
        % plot learning results for one trajectory
        function plot_learning(obj,trj)
            
            trj.plot_learning_output(obj);
            trj.plot_learning_input(obj);
            
        end
    end
end