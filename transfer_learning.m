%% Transfer Learning 
% from one trajectory to another
% comparing TGP with MPC, ILC and TGP from scratch
%
% AUTHOR: Okan Koc, ETHZ IDSC Lab, Zurich, October 2012 - January 2014
%
% -------------------------------------------------------------------------

%% Prepare experiments and estimate hyperparameters     
clc; clear; close all;
fix_root();
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

%% Generate first trajectory along which TGP runs
shape = 'Wave';
incr_y = [0.126, 0.145, 0.309];
incr_z = [0.1833, 0.412, 0.4691];
y_coord = ([0 1.2 2.4 3.2] + [0 incr_y])/10;
z_coord = ([0 1.8 0.0 1.2] + [0 incr_z])/10;
spline = [y_coord; z_coord];
trj1 = quad.trajectory(shape, spline);

%% Generate second trajectory along which TGP runs
shape = 'Wave';
% determine randomized coordinates of the new wave
incr_y = 0.5*rand(1,3);
incr_z = 0.5*rand(1,3);
%{
%rotate wave by alpha degrees
alpha = 0;
rot = [cos(alpha), sin(alpha);
       -sin(alpha), cos(alpha)];
newWave = rot * [y_coord; z_coord];
y_coord = newWave(1,:);
z_coord = newWave(2,:);
%}
y_coord = ([0 1.2 2.4 3.2] + [0 incr_y])/10;
z_coord = ([0 1.8 0.0 1.2] + [0 incr_z])/10;
spline = [y_coord; z_coord];
trj2 = quad.trajectory(shape, spline);

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
num_epi = 1;

%% Run TGP on first trajectory
% initialize TGP
tgp = TGP(quad,est,false);

for iter = 1:num_epi
    % start on the trajectory
    dev_gp = 0; i = 0;
    x_gp = trj1.s(:,1); 
    u_gp = trj1.unom(:,1);
    while i < trj1.N-1 %&& dev_gp < tgp.epsilon
        i = i+1;
        % apply control signal
        u_gp(:,i) = tgp.control(i,trj1,quad,x_gp(:,i));
        % observe output by evolving system
        x_gp(:,i+1) = quad.evolve(trj1.t(i),x_gp(:,i),u_gp(:,i));
        % get deviation for tgp
        dev_gp = quad.COST.fnc(x_gp(:,i+1),trj1.s(:,i+1));
        tgp.update_data(quad,dev_gp,u_gp(:,i));
        % display every 20th iteration if process is slow
        if(rem(i,20) == 0), fprintf('Time stage %d ...\n', i); end
    end
    % display stopping time
    fprintf('Iteration %d stopping time %d \n',iter,i);
    % record the performance
    trj1.addPerformance(u_gp,x_gp,[],quad.COST,tgp);        
end

%% Run TGP on second trajectory
% initialize TGP from scratch
tgp_scratch = TGP(quad,est,false);

for iter = 1:num_epi
    % start on the trajectory
    [dev_gp1,dev_gp2] = deal(0); i = 0;
    [x_gp1,x_gp2] = deal(trj2.s(:,1));
    [u_gp1,u_gp2] = deal(trj2.unom(:,1));
    while i < trj2.N-1 %&& dev_gp < tgp.epsilon
        i = i+1;
        % apply control signal
        u_gp1(:,i) = tgp.control(i,trj2,quad,x_gp1(:,i));
        u_gp2(:,i) = tgp_scratch.control(i,trj2,quad,x_gp2(:,i));
        % observe output by evolving system
        x_gp1(:,i+1) = quad.evolve(trj2.t(i),x_gp1(:,i),u_gp1(:,i));
        x_gp2(:,i+1) = quad.evolve(trj2.t(i),x_gp2(:,i),u_gp2(:,i));
        % get deviation for tgp
        dev_gp1 = quad.COST.fnc(x_gp1(:,i+1),trj2.s(:,i+1));
        dev_gp2 = quad.COST.fnc(x_gp2(:,i+1),trj2.s(:,i+1));
        tgp.update_data(quad,dev_gp1,u_gp1(:,i));
        tgp_scratch.update_data(quad,dev_gp2,u_gp2(:,i));
        % display every 20th iteration if process is slow
        if(rem(i,20) == 0), fprintf('Time stage %d ...\n', i); end
    end
    % display stopping time
    fprintf('Iteration %d stopping time %d \n',iter,i);
    % record the performance
    trj2.addPerformance(u_gp1,x_gp1,[],quad.COST,tgp);        
    trj2.addPerformance(u_gp2,x_gp2,[],quad.COST,tgp_scratch);        
end

%% Run ILC on second trajectory
x_ilc = trj2.PERF(1).x;
ilc = ILC(quad,trj2);
dev_ilc = zeros(1,trj2.N);
for i = 1:num_epi
    % generate feedforward controls
    u_ilc = ilc.control(trj2,quad,x_ilc - trj2.s);
    x_ilc = quad.evolve_full(trj2.t,x_ilc,u_ilc);
    trj2.addPerformance(u_ilc,x_ilc,[],quad.COST,ilc);
end    

%% Plot the performances

costs = [ilc.sse,tgp_scratch.sse,tgp.sse(2)];
figure;
bar(costs);