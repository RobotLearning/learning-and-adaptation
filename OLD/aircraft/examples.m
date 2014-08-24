%% EXAMPLES FOR LEARNING HYPERPARAMETERS
%--------------------------------------------------------------------------
% Called once by the GP main script (GP_Main)
%
% Author: Okan Koc, IDSC Lab, ETH Zurich
%
% -------------------------------------------------------------------------

% Run n-many different trajectories
% Supervise with n-many nominal control inputs and costs
% Possible shapes are n different sinusoidals
n = 5;
% initialize variables
N = 100; % number of time stages in each example
% cells are used because N might not be the same in all models
us = cell(n,1); 
x_noms = cell(n,1);
x_reals = cell(n,1);

% constructing z1 and z2
r = 10*ones(1,n);
w = (1:n)*(pi/10);
w(n) = 2*pi; % last try is blowing up!
t = linspace(h,h*N,N);

for i = 1:n
    z1 = r(i) - r(i)*cos(w(i)*t); 
    z2 = r(i)*sin(w(i)*t);
    ddz1 = (w(i)^2)*r(i)*cos(w(i)*t);
    ddz2 = -(w(i)^2)*r(i)*sin(w(i)*t);
    traj = [z1; z2; ddz1; ddz2];
    [t,u_trj,x_nom] = traj_gen(traj,PAR,h);
    us{i} = u_trj;
    x_noms{i} = x_nom;
    % plot nominal trj and deviation
    x_real = sim_RK4(t,x_nom(:,1),us{i},0,PAR,fun_real);
    plot_nominal_real_trj;
    % Covariance matrix for sensor error / noise model
    Omega = eps * eye(N*dim_x);
    
    %%%%%%%%%%%%%% ADD NOISE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    % vectorize x 
    x_vec = x_real(:);
    % add noise w with covar Omega
    x_vec_noisy = x_vec + chol(Omega)*randn(length(x_vec),1);
    % back to matrix form
    x_noisy = reshape(x_vec_noisy,dim_x,N);
    x_reals{i} = x_noisy; %IC is also noisy
    %}
    x_reals{i} = x_real;
end
close all;

%%%%%%%%%%%%%%%%%%%%%%%% GET CONTEXTS AND COSTS %%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize costs and contexts cells
costs = cell(n,1);
contexts = cell(n,1);
% prepare costs
for i = 1:n
    %initialize contexts and costs
    context = zeros(dim_ctx, N-1);
    cost = zeros(N-1, 1);
    for j = 1:N-1
        % as context, we keep only present state + future state
        context(:,j) = [x_reals{i}(:,j); x_noms{i}(:,j+1)];
        % get the 2-norm difference squared
        dev = x_reals{i}(:,j+1) - x_noms{i}(:,j+1);
        cost(j) = dev'*Q*dev;
        % subtract nominal model prediction error
        next = step_RK4(h,x_reals{i}(:,j),us{i}(:,j),CON,PAR,fun_nom);
        % get nominal deviation
        pred_dev = next - x_noms{i}(:,j+1);
        pred_err = pred_dev'*Q*pred_dev;
        cost(j) = cost(j) - pred_err;
    end
    contexts{i} = context;
    costs{i} = cost;
end

% number of total training samples
num_tot = (N-1)*n;