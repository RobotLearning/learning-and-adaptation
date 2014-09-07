%% SET PARAMETERS
%--------------------------------------------------------------------------
% Called once by the GP main script (GP_Main)
%
% Author: Okan Koc, IDSC Lab, ETH Zurich
%
% -------------------------------------------------------------------------

%%%%%%%%%%%%% Define System and Constraint Structures %%%%%%%%%%%%%%%%%%%%%

% System constants
PAR.r = 100.0; % m
PAR.m = 1000; % kg
PAR.J = 100; % kg m2
PAR.g = 9.81; % m/s2
PAR.Craft.A = 100; %m2
% Wind disturbance
PAR.Wind.Pressure = 100; % N/m2
PAR.Wind.Angle = 0 * pi/180; % wind angle
% No constraints
CON = Inf;
CON.F1 = 1e04;
CON.F2 = 2e04;
% bounds for the convex space in CGP-UCB
r = 1000;
% flags for the input bounds/derivative bounds
flag_bnd = false;
flag_drv_bnd = false;

%%%%%%%%%%%%%%%%%%%%%%% Simulation Values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dimension of the x vector
dim_x = 6;
% dimension of the action space
dim = 2;
% only some dimensions are tracked
dims_tr = [1,3]; % x and y 
% time step h 
h = 0.02;
% noise and initial error
eps = 0.00003;
% set horizon size for calculating converging traj - CGP horizon, not MPC
horizon = 5;
% Scaling matrix Sw
Sw = eye(dim_x);
%Sw = diag([1, 1, 1, 1, 0.6, 0.5]); 
Q = Sw'*Sw;
% dimension of the context space : current state + desired state
dim_ctx = dim_x * 2; 
% real function handle for calculating 'real' trajectories
fun_real = @aircraftRealDynamics;
% nominal model to generate nominal prediction costs
fun_nom = @aircraftNominalDynamics;