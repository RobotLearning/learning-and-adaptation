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

% desired amplitudes and frequencies of the sinusoidal trajectories
A_d = linspace(1,5,n);
w_d = 2*ones(1,n);
time = linspace(h,1,N);

for i = 1:n
    trj = A_d(i) * sin(w_d(i)*2*pi*time);
    dtrj = w_d(i)*2*pi * A_d(i) * cos(w_d(i)*2*pi*time);
    traj = [trj; dtrj];
    [t,u_trj,x_nom] = traj_gen(traj,PAR,h);
    us{i} = u_trj;
    x_noms{i} = x_nom;
    % plot nominal trj and deviation
    x_real = sim_RK4(t,x_nom(:,1),us{i},0,PAR,fun_real);
    plot_nominal_real_trj;
    % Covariance matrix for sensor error / noise model
    Omega = eps * eye(N*dim_x);
    
    %%%%%%%%%%%%%% ADD NOISE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % vectorize x 
    %x_vec = x_real(:);
    % add noise w with covar Omega
    %x_vec_noisy = x_vec + chol(Omega)*randn(length(x_vec),1);
    % back to matrix form
    %x_noisy = reshape(x_vec_noisy,dim_x,N(i));
    %x_reals{i} = x_noisy; %IC is also noisy
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