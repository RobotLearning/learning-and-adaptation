% Class implementing a Gaussian Process regression model

classdef gp < handle
    
    properties          

        % Hyperparameters
        hp
        % Training data
        x,y
        % Covariance matrix between (training) data points
        cov
        % Kernel function (handle)
        kernel
        
    end
    
    methods (Access = public)
        
        %% Initialize data points, kernel structure and covariance matrix
        function obj = gp(hp,xs,ys)
            
            obj.x = xs;
            obj.y = ys;
            obj.hp = hp;
            obj.set_kernel(hp);
            obj.build_cov(xs);
            
        end
        
        % Estimate hyperparameters using maximum likelihood
        function fit_hp(obj,hp0)
             
            try
                addpath(genpath('../gpml-matlab-v3.1-2010-09-27'));
            catch
                error('GPML toolbox not found!');
            end
            
            covar = @covSEard;
            lik = @likGauss;
            inf = @infExact;
            mean = [];
            cov_hp_len = size(obj.x,1) + 1; % scale hp added
            hyp.mean = [];
            hyp.cov = [log(hp0.l);log(sqrt(hp0.scale))];
            hyp.lik = log(sqrt(hp0.noise.var));            
            Ncg = 1000; % 100 line search steps
            hyp = minimize(hyp, @gp, Ncg, ...
                              inf, mean, covar, lik, obj.x', obj.y);
                          
            % update hp field after hp estimation
            obj.hp.l = exp(hyp.cov(1:end-1));
            obj.hp.scale = exp(hyp.cov(end))^2;
            obj.hp.noise.var = exp(hyp.lik)^2;
            obj.set_kernel(obj.hp);
        end
        
        
        %% predict/update mean and variance at test point(s)
        function [mu,s2] = predict(obj,xstar)
            
            % calculate kernel vector k and store
            vec = obj.cov_test_and_data(xstar);
            % create mu function 
            mu = vec' * (obj.cov \ obj.y);
            % initial covar is 1
            kxx = obj.kernel(xstar,xstar);
            % subtract information gain
            s2 = kxx - vec' * (obj.cov \ vec);
        end
        
        function [mu,s2] = predict_mesh(obj,mesh)
           
            % calculate kernel matrix K at mesh test points
            cov_mesh = obj.cov_mesh(mesh);
            m = size(mesh,2);
            n = length(obj.y);
            if n == 0
                mu = zeros(1,m);
                s2 = cov_mesh;
            else
                mat = zeros(n,m);
                for i = 1:m
                    vec = obj.cov_test_and_data(mesh(:,i));
                    mat(:,i) = vec;
                end
                mu = mat' * (obj.cov \ obj.y);
                s2 = cov_mesh - mat' * (obj.cov \ mat);
            end
            s2 = (s2 + s2')/2; % to enforce symmetric matrix
        end
        
        % update mean and variance
        function update(obj,xnew,ynew)
            
            obj.update_cov(xnew);
            obj.x = [obj.x, xnew];
            obj.y = [obj.y; ynew];
        end        
    end
    
    methods (Access = private)
        
        %% Kernel related methods
        
        % Kernel called by the ker_matrix, ker_matrix_iter and ker_vector
        % functions.
        % x1 - Vector 1
        % x2 - Vector 2 
        %
        % Sets up a Gaussian or linear kernel function depending on 
        % the type field and hyperparameters
        %
        % IMPORTANT: does not consider noise as part of kernel
        %
        function set_kernel(obj,hp)

            L = hp.l;
            type = hp.type;
            s = hp.scale;
            switch type
                case 'squared exponential iso'
                    if length(L) == 1
                        obj.kernel = @(x1,x2) s * exp(-(norm(x1(:)-x2(:),2)^2)/(2*(L^2)));
                    else
                        error('Lengthscale parameter should be scalar!');
                    end
                case 'linear iso'
                    if isempty(L)
                        obj.kernel = @(x1,x2) s * x1(:)'*x2(:);
                    else
                        error('Lengthscale parameter should be empty!');
                    end
                case 'squared exponential ard'
                    InvGamma = diag(1./(L.^2));
                    obj.kernel = @(x1,x2) s * exp(-0.5*((x1(:)-x2(:))')*InvGamma*(x1(:)-x2(:)));
                case 'linear ard'
                    InvGamma = diag(1./(L.^2));
                    obj.kernel = @(x1,x2) s*x1(:)'*InvGamma*x2(:);
                otherwise
                    error('Unrecognized kernel type');
            end
        
        end        
        
        % Constructs the kernel matrix all at once.
        % Only used to generate test functions.
        % INPUTS: 
        % x - training data points
        % kernel - covariance function
        % OUTPUTS
        % out - matrix of size nxn, where n is the length of us
        function build_cov(obj,x)

            len = size(x,2);
            Kmat = zeros(len,len);
            for i = 1:len
                for j = i:len
                    Kmat(i,j) = obj.kernel(x(:,i), x(:,j));
                end
            end
            obj.cov = Kmat + Kmat' - diag(diag(Kmat));
            obj.cov = obj.cov + obj.hp.noise.var * eye(len);
        end
        
        % Construct the kernel matrix iteratively
        % INPUTS: 
        % t - last time stage (up to horizon)
        % us - control law points 
        % fun - function handle for the kernel
        % mat - previous kernel matrix
        % OUTPUT
        % Kmat - updated kernel matrix at the previous points and current test pt.
        function update_cov(obj,xstar)

            % updates iteratively
            k = obj.cov_test_and_data(xstar);
            Kmat = obj.cov;
            % since kernel matrix is symmetric 
            % simply update the last row with 
            % the last column
            Kmat(:,end+1) = k;
            kxx = obj.kernel(xstar,xstar) + obj.hp.noise.var;
            Kmat(end+1,:) = [k;kxx]';
            obj.cov = Kmat;
        end
        
        % construct the kernel vector - covariances of the test pt and previous pts
        % INPUTS: 
        % xstar - test point
        % x     - data points 
        % OUTPUT
        % cov - kernel column vector describing the covariances between
        %       the previous points and current test pt.
        %
        function cov = cov_test_and_data(obj,xstar)

            n = size(obj.x,2);
            cov = zeros(n,1);
            for i = 1:n
                cov(i) = obj.kernel(obj.x(:,i),xstar);
            end
        
        end
        
        % Construct the kernel matrix between the 
        % mesh of test points
        %
        function cov = cov_mesh(obj,xs)
           
            len = size(xs,2);
            cov = zeros(len);
            for i = 1:len
                for j = i:len
                    cov(i,j) = obj.kernel(xs(:,i), xs(:,j));
                end
            end
            cov = cov + cov' - diag(diag(cov));
        end
        
        
    end

end