%% EXAMPLES FOR LEARNING HYPERPARAMETERS
%--------------------------------------------------------------------------
% Called once by the GP main script (GP_Main)
%
% Author: Okan Koc, IDSC Lab, ETH Zurich
%
% -------------------------------------------------------------------------

% Run n-many different trajectories
% Supervise with n-many nominal control inputs and costs
% Possible shapes are as follows: 
% shapes = {'S shape', 'Horizontal', 'Diagonal', 'Arc', 'Wave', 'Circle', ...
%           'Question mark', 'A shape'};
shapes = {'Wave','Wave','Wave','Wave','Wave'};
n = length(shapes);
% determine coordinates of the wave
y_coord = [0 1.0 2.0 3.0;
           0 1.2 1.8 2.5;
           0 0.8 2.2 2.7;
           0 0.6 2.4 3.0;
           0 1.0 2.5 3.5];
% to reduce trajectory time
y_coord = y_coord/10;
z_coord = [0 1.5 0.0 1.5;
           0 1.2 0.2 1.5;
           0 1.8 0.4 1.8;
           0 0.5 1.2 1.8;
           0 2.0 1.0 0.0;];
z_coord = z_coord/10;

% dont run code all the time
%matfile = 'ctx_quad_actuator.mat';
matfile = 'ctx_quad_wind.mat';
savefile = strcat(qcopterFolder, '/', matfile);

if ~exist(savefile,'file')

% initialize variables
N = zeros(n,1); % number of time stages in each example
us = cell(n,1);
x_noms = cell(n,1);
x_reals = cell(n,1);

for i = 1:n
    traj = shapes{i};
    [t,u_trj,x_nom] = traj_gen(traj,PAR,CON,h,y_coord(i,:),z_coord(i,:),i);
    us{i} = u_trj;
    x_noms{i} = x_nom;
    % Covariance matrix for sensor error / noise model
    N(i) = length(t);
    Omega = eps * eye(N(i)*dim_x);
    x_real = sim_RK4(t,x_nom(:,1),us{i},CON,PAR,fun_real);
    % plot nominal trj and deviation
    plot_nominal_real_trj;

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
    context = zeros(dim_ctx, N(i)-1);
    cost = zeros(N(i)-1, 1);
    for j = 1:N(i)-1
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

%%{
save(savefile', 'us', 'contexts', 'costs', 'N');
else
load(savefile);
end
%}

% number of total training samples
num_tot = sum(N) - n;