%% SET PARAMETERS
%--------------------------------------------------------------------------
% Called once by the GP main script
%
% Author: Okan Koc, IDSC Lab, ETH Zurich
%
% -------------------------------------------------------------------------

%%%%%%%%%%%%% Define System and Constraint Structures %%%%%%%%%%%%%%%%%%%%%

% System constants
PAR.k = 1.0;
PAR.m = 1.0;
PAR.c = 1.0;
% state-space representation of the nominal model
PAR.A = [0, 1; -PAR.k/PAR.m, -PAR.c/PAR.m]; 
PAR.B = [0; 1/PAR.m];
% No constraints
CON = Inf;
% bounds for the convex space in CGP-UCB
r = 100;
% flags for the input bounds/derivative bounds
flag_bnd = false;
flag_drv_bnd = false;

%%%%%%%%%%%%%%%%%%%%%%% Simulation Values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dimension of the x vector
dim_x = 2;
% dimension of the action space
dim = 1;
% all dimensions are tracked
dims_tr = 2; 
% time step h 
h = 0.01;
% noise and initial error
eps = 0.00003;
% set horizon size for calculating converging traj - CGP horizon, not MPC
horizon = 5;
% Scaling matrix Sw
%Sw = eye(dim_x);
Sw = diag([1, 1]); 
Q = Sw'*Sw;
% dimension of the context space : current state + desired state
dim_ctx = dim_x * 2; 
% real function handle for calculating 'real' trajectories
fun_real = @springDamperRealDynamics;
% nominal model to generate nominal prediction costs
fun_nom = @springDamperNominalDynamics;