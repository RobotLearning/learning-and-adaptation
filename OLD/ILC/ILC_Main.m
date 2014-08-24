%% Set system parameters and constraints

clc; clear, close all;

%%%%%%%%%%%%% Define System and Constraint Parameters %%%%%%%%%%%%%%%%%%%%%

%System constants
PAR.Iy = 0.0023;
PAR.L = 0.17;
PAR.m = 0.468;
PAR.g = 9.81;
%Disturbance constants
PAR.Quad.A = 0.4 * 0.4; % 40 cm by 40 cm quadrocopter
PAR.Fan.Angle = 0 * pi/180; % 45 degrees blowing northwest
PAR.Fan.Length = 0.4; % 40 cm fan
PAR.Fan.GroundPos = 1.3; % Fan centre at 1.3 m in y direction
PAR.Fan.BlowPressure = 50; % Fan blow pressure N/m2
%Single thrust constraints
CON.fi_max = 4.5;
CON.fi_min = 0.4;
CON.fi_dot_max = 27;
%Collective thrust constraints
CON.fmax = 4*CON.fi_max;
CON.fmin = 4*CON.fi_min;
%Angle constraints
CON.phi_max = pi/2;
CON.phi_dot_max = 22;
CON.phi_ddot_max = 150;
%Number of constraint evaluation points
CON.NumT = 500;

%%%%%%%%%%%%%%%%%%%%%%% Simulation Values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dimension of the x vector
dim_x = 5;
% dimension of control input
dim = 2;
% time step h 
h = 0.02;
% noise and initial error
eps = 0.003;
eps_M = 0.05; % covar of the lifted-domain noise
d0_hat = 0;

%% Create desired trajectory and linearized model

% generate trajectory from script
traj = 'Wave';
[t,u_trj,x,fc,rollrate] = traj_gen(traj,PAR,CON,h);

%%%%%%%%%%%%% INCREASE CONSTRAINTS for learning process %%%%%%%%%%%%%%%%%%%
%Single thrust constraints
CON.fi_max = 5.5;
CON.fi_min = 0.25;
CON.fi_dot_max = 51;
%Collective thrust constraints
CON.fmax = 4*CON.fi_max;
CON.fmin = 4*CON.fi_min;
%Angle constraints
CON.phi_dot_max = 25;
CON.phi_ddot_max = 200;

%Add to plots increased constraints
%{
figure(2);
subplot(2,1,1)
fmax = plot([t(1),t(end)],[CON.fmax,CON.fmax],'--r');
legend([fc,fmax],'Thrust','Constraints')

subplot(2,1,2)
phi_dot_max = plot([t(1),t(end)],[CON.phi_dot_max,CON.phi_dot_max],'--r');
legend([rollrate,phi_dot_max], 'Roll rate', 'Constraints');
%}

%%%%%%%%%%%%%%%% LINEARIZED MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct linearized model around the desired trajectory
% structure containing ILC matrices returned
STR = construct_model(t,x,PAR,u_trj,CON);

% initial Kalman filter variance
N = length(t);
Nu = N-1; % u input needs one less
P0 = eye(dim_x*N);
% lifted-domain noise covariance
M = eps_M * eye(N*dim_x); % covariance of process x measurement noise
Omega0 = eps * eye(N*dim_x); % initial covariance of d noise

%% Iterative Learning Control

w = ones(1,N+1); % equal weighting
% SSE-cost for trajectory deviation
sse_cost = @(x1,x2,w) sum(sum((x1([1,3],:)-x2([1,3],:)).^2));

%%%%%%%%%%%%%%%%% CONSTRUCT SCALING MATRIX S %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% penalize only deviations given by the matrix
Sx = diag([2, 1, 2, 1, 5]); % scaling matrix
Sw = diag([1, 0.1, 1, 0.1, 0]); % weighing trajectory deviation by these values
%Sx = eye(dim_x);
%Sw = diag([1, 0, 1, 0, 0]);
w_l = 0.0; % put more/less weight on last w_l second 
% emphasizing the performance objective of reaching upright position
last = w_l * round(1/h);
Sx1 = cell(1,N - last);
[Sx1{:}] = deal(Sx * Sw);
Sx1 = blkdiag(Sx1{:});
Sx2 = cell(1,last);
% termination matrix
Tw = diag(ones(1,dim_x));
mat = Tw * Sx * Sw;     
[Sx2{:}] = deal(mat);
Sx2 = blkdiag(Sx2{:});
S = blkdiag(Sx1,Sx2);

%%%%%%%%%%%%%%%% PASS TO STRUCTURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pass necessary parameters through structure
STR.dim_x = dim_x;
STR.dim = dim;
STR.COV.Omega = Omega0;
STR.COV.M = M;
STR.CON = CON;
STR.PAR = PAR;
STR.d0 = d0_hat * ones(N*dim_x,1);
STR.P0 = P0;
STR.w = w;
STR.sse_cost = sse_cost;
STR.S = S;

%%%%%%%%%%%%% ITERATIVE LEARNING CONTROL LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fun = @quadrocopterError;
[err, u_app] = ILC(t,x,u_trj,fun,STR,5);
% errors can be reused in some other code