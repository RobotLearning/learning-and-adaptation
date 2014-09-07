% Trajectory Tracking using Gaussian Processes (TGP)
% Comparing TGP with ILC and MPC
%
% AUTHOR: Okan Koc, ETHZ IDSC Lab, Zurich, October 2012 - January 2014
%
% -------------------------------------------------------------------------

%% Prepare experiments and estimate hyperparameters     
clc; clear; close all;
fix_root('~/');
quad_initialize; % create structures for quadrotor operation here
quad = Quadrotor(par,dist,con,Q);
shapes = {'Wave','Wave','Wave','Wave','Wave'};
coord = {[0 1.0 2.0 3.0; 0 1.5 0.0 1.5]/10, ...
         [0 1.2 1.8 2.5; 0 1.2 0.2 1.5]/10, ...
         [0 0.8 2.2 2.7; 0 1.8 0.4 1.8]/10, ...
         [0 0.6 2.4 3.0; 0 0.5 1.2 1.8]/10, ...
         [0 1.0 2.5 3.5; 0 2.0 1.0 0.0]/10};
quad.flg_exp = false;
quad.flg_est = false;     
est = quad.experiment(shapes,coord);

%% Generate trajectory along which TGP runs
shape = 'Wave';
incr_y = [0.126, 0.145, 0.309];
incr_z = [0.1833, 0.412, 0.4691];
y_coord = ([0 1.2 2.4 3.2] + [0 incr_y])/10;
z_coord = ([0 1.8 0.0 1.2] + [0 incr_z])/10;
spline = [y_coord; z_coord];
trj = quad.trajectory(shape, spline);

%% Increase constraints for the learning process
%Single thrust constraints
CON.fi_max = 5.5;
CON.fi_min = 0.25;
CON.fi_dot_max = 51;
%Collective thrust constraints
CON.fmax = 4*CON.fi_max;
CON.fmin = 4*CON.fi_min;
CON.f_dot_max = 4*CON.fi_dot_max;
%Angle constraints
CON.phi_dot_max = 25;
CON.phi_ddot_max = 200;
%These don't change
CON.phi_max = pi/2; 
CON.NumT = 500;
quad.CON = CON;

% update bounds
mat = [CON.fmin, CON.fmax; 
      -CON.phi_dot_max, CON.phi_dot_max];
quad.bound = mat;

% total number of episodes
num_epi = 6;

%% Learning process for TGP
% initialize tgp
tgp = TGP(quad,est,true);
for iter = 1:num_epi
    % start on the trajectory
    dev_gp = 0; i = 0;
    x_gp = trj.s(:,1); 
    u_gp = trj.unom(:,1);
    while i < trj.N-1 && dev_gp < tgp.epsilon
        i = i+1;
        % apply control signal
        u_gp(:,i) = tgp.control(i,trj,quad,x_gp(:,i));
        % observe output by evolving system
        x_gp(:,i+1) = quad.evolve(trj.t(i),x_gp(:,i),u_gp(:,i));
        % get deviation for tgp
        dev_gp = quad.COST.fnc(x_gp(:,i+1),trj.s(:,i+1));
        tgp.update_data(quad,dev_gp,u_gp(:,i));
        % display every 20th iteration if process is slow
        if(rem(i,20) == 0), fprintf('Time stage %d ...\n', i); end
    end
    % display stopping time
    fprintf('Iteration %d stopping time %d \n',iter,i);
    % record the performance
    trj.addPerformance(u_gp,x_gp,[],quad.COST,tgp);        
end

%% Run ILC
x_ilc = trj.PERF(1).x;
ilc = ILC(quad,trj);
dev_ilc = zeros(1,trj.N);
for i = 1:num_epi
    % generate feedforward controls
    u_ilc = ilc.control(trj,quad,x_ilc - trj.s);
    x_ilc = quad.evolve_full(trj.t,x_ilc,u_ilc);
    trj.addPerformance(u_ilc,x_ilc,[],quad.COST,ilc);
end    

%% Model Predictive Control
% initialize mpc
horizon = 5; mpc = MPC(quad,horizon);
for iter = 1:num_epi
    % start on the trajectory
    x_mpc = trj.s(:,1); 
    u_mpc = trj.unom(:,1);
    for i = 1:trj.N - 1
        % apply control signal 
        u_mpc(:,i) = mpc.control(i,trj,quad,x_mpc(:,i));    
        % observe output by evolving system
        x_mpc(:,i+1) = quad.evolve(trj.t(i),x_mpc(:,i),u_mpc(:,i)); 
        % display every 20th iteration if process is slow
        if(rem(i,20) == 0), fprintf('Time stage %d ...\n', i); end
    end
    % record the performance
    trj.addPerformance(u_mpc,x_mpc,[],quad.COST,mpc);
end

%% Plot the performances
quad.plot_learning(trj);
save('trj.mat','trj');

% plot the costs
figure;
plot(1:ilc.epi,ilc.sse,ilc.color,...
    1:mpc.epi,mpc.sse,mpc.color,...
    1:tgp.epi,tgp.sse,tgp.color);