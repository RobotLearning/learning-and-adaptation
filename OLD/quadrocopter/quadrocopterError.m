function [traj, err] = quadrocopterError(t, u_trj, u, x, PAR, COV, CON)

% calculates the trajectory and error given the dynamics

dim_x = size(x,1);
% current control input u
u_cur = u_trj + u;
% state trajectories for this iteration
traj = x;
Omega = COV.Omega(dim_x+1:end,dim_x+1:end);
M = COV.M(dim_x+1:end,dim_x+1:end);

% simulate real trajectory with RK4
handle = @quadrocopterRealDynamics;
traj = sim_RK4(t,traj(:,1),u_cur,CON,PAR,handle);

% vectorize x_iter into N*dim dimensions
x_iter_vec = traj(:);
x_iter_vec = x_iter_vec(dim_x+1:end); % start from x(2)
% add lifted-domain error noise w with covariance Omega
x_vec_noisy = x_iter_vec + chol(Omega)*randn(length(x_iter_vec),1);
% add observation noise with covariance M
y = x_vec_noisy + chol(M)*randn(length(x_iter_vec),1);

% find deviation in observation
x_vec = x(:); 
x_vec = x_vec(dim_x+1:end);
err = y - x_vec;
% start from x(1)
err = [zeros(dim_x,1); err];