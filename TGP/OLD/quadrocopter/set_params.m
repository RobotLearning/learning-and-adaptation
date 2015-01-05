%% SET PARAMETERS
%--------------------------------------------------------------------------
% Called once by the GP main script (GP_Main)
%
% Author: Okan Koc, IDSC Lab, ETH Zurich
%
% -------------------------------------------------------------------------

%%%%%%%%%%%%% Define System and Constraint Structures %%%%%%%%%%%%%%%%%%%%%

%System constants
PAR.Iy = 0.0023;
PAR.L = 0.17;
PAR.m = 0.468;
PAR.g = 9.81;
%Disturbance constants
PAR.Quad.A = 0.4 * 0.4; % 40 cm by 40 cm quadrocopter
PAR.Fan.Angle = 0; % angle of the blowing disturbance
PAR.Fan.Length = 0.4; % 40 cm fan
PAR.Fan.GroundPos = 1.3; % Fan centre at 1.3 m in y direction
PAR.Fan.BlowPressure = 60; % Fan blow pressure N/m2
%Single thrust constraints
CON.fi_max = 4.5;
CON.fi_min = 0.4;
CON.fi_dot_max = 27;
%Collective thrust constraints
CON.fmax = 4*CON.fi_max;
CON.fmin = 4*CON.fi_min;
%Angle constraints
CON.phi_max = pi/2; % only used for trajectory generation
CON.phi_dot_max = 22;
CON.phi_ddot_max = 150;
%Number of constraint evaluation points
CON.NumT = 500;
% bounds of the cube
r = max([CON.fmax; CON.phi_dot_max]); 
% flags for the input bounds/derivative bounds
flag_bnd = true;
flag_drv_bnd = false;

%%%%%%%%%%%%%%%%%%%%%%% Simulation Values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dimension of the x vector
dim_x = 5;
% dimension of the action space
dim = 2;
% 2 dimensions, y and z, are tracked only
dims_tr = [1,3]; 
% time step h 
h = 0.02;
% noise and initial error
eps = 0.00003;
% set horizon size for calculating converging traj - CGP horizon, not MPC
horizon = 5;
% Scaling matrix Sw
%Sw = eye(dim_x);
Sw = diag([1, 1, 1, 1, 0]); 
Q = Sw'*Sw;
% dimension of the context space : current state + desired state
dim_ctx = dim_x * 2; 
% real function handle for calculating 'real' trajectories
fun_real = @quadrocopterRealDynamics;
% nominal model to generate nominal prediction costs
fun_nom = @quadrocopterNominalDynamics;