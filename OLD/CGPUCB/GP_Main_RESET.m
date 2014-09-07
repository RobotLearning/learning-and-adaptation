%% CGP-UCB script for a dynamical system simulation
%
% The method resets here if the deviation is above a certain threshold.
%
% AUTHOR: Okan Koc, ETHZ IDSC Lab, Zurich, October 2012 - March 2013
%
% -------------------------------------------------------------------------

clc; clear; close all;

% set the shared folders
root = [cd, '\']; % in windows
%root = '/home/okan/Learning-and-Adaptation/';
cgpFolder = [root, 'CGPUCB'];
mpcFolder = [root, 'MPC'];
% set dynamics here
qcopterFolder = [root, 'quadrocopter'];
springFolder = [root, 'springDamper'];
aircraftFolder = [root, 'aircraft'];
pendulumCartFolder = [root, 'pendulumOnCart'];
pendulumSwingUpFolder = [root, 'pendulumSwingUP'];
robotKinematicsFolder = [root, 'robotTwoWheelsKinematics'];
% remove unnecessary model folder paths
rmpath(springFolder); rmpath(aircraftFolder); rmpath(pendulumCartFolder);
rmpath(pendulumSwingUpFolder); rmpath(robotKinematicsFolder);
rmpath(qcopterFolder);
% add the folder you want
addpath(cgpFolder);
addpath(mpcFolder);
addpath(qcopterFolder);

%% Set system parameters
set_params;

%% Generate examples, contexts and costs
examples;

%% Fit hyperparameters on nominal model prediction error
fit_hyperparameters;
%debug;

%% Generate first trajectory along which GP-UCB runs
generate_trj;

%% Prepare variables for GP-UCB
prepare_gpucb;

%%%% FOR DEBUGGING 
%STR.handle = fun_real;
% each run learns from the previous runs
runs = 3;
STR.MPC.HORIZON = 1;
STR.GPMPC.HORIZON = 5;
%STR.CGP.HORIZON = 1; % used to find converging trj.

%% Bandit process

%%%%%%%%%%%%%%%%% DEFINE RESET THRESHOLD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshold = 0.03;

%%%%%%%%%%%%%%%%% LEARN THROUGH MULTIPLE COMPLETE TRIALS %%%%%%%%%%%%%%%%%%

% INVESTIGATING OVERFITTING EFFECTS
past_example_data = num_tot;
window = 0;
flag_GPMPC = false;
% initialize kernel for the gp-ucb-original
STR.Kinv = 1/(1+exp(GPSTR.hyp.lik)^2);
STR.K = 0;
data = zeros(size(x_nom,1), size(x_nom,2), runs);

for run = 1:runs
    
    % adjust x by adding endpt to the horizon
    x_hor = [x_nom, repmat(x_nom(:,end),1,horizon)];    
    % start from the initial point
    x_GP(:,1) = x_nom(:,1); 
    x_MPC(:,1) = x_nom(:,1);
    u(:,1) = u_trj(:,1);
    u_mpc = zeros(dim,length(x_nom));
    
    fprintf('Run number %d ...\n', run);
            
    %%%%%%%%%%%%%%%%% LEARN THROUGH ONE TRIAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    for t = 1:20
        cgp_ucb;         
        call_mpc;
        
        dev = norm(x_GP(dims_tr,t+1) - x_nom(dims_tr,t+1));
        % break if trajectory deviation exceeds threshold
        if dev >= threshold, break; end

        % display every 10th iteration if process is slow
        if(rem(t,10) == 0), fprintf('Time stage %d ...\n', t); end
    end
    
    % display runtime
    fprintf('Trajectory simulation took %f seconds. \n', toc);
    
    % show results
    show_results;
    
    %%%%%%%%%%%%%%%%%%%% CHANGE TRAJECTORY HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if run < runs
        change_trj;
    end
    
    % save result 
    %data(:,:,run) = x_GP;
end

%{
% plot results for presentation
figure;
plot(nom1(dims_tr(1),:), nom1(dims_tr(2),:), '-r');
hold on;
for i = 1:runs-1
    p = plot(data(dims_tr(1),:,i), data(dims_tr(2),:,i));
    set(p,'Color',rand(1,3));
end
legend('Nom', 'It 1', 'It 2');
%}