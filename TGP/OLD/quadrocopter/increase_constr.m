%%%%%%%%%%%%% INCREASE CONSTRAINTS FOR LEARNING PROCESS %%%%%%%%%%%%%%%%%%%

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

% Establish bounds for variables
STR.MPC.bounds = [CON.fmin, CON.fmax; -CON.phi_dot_max, CON.phi_dot_max];
STR.CGP.bounds = [CON.fmin, CON.fmax; -CON.phi_dot_max, CON.phi_dot_max];
STR.GPMPC.bounds = [CON.fmin, CON.fmax; -CON.phi_dot_max, CON.phi_dot_max];

% needed to calculate contexts (opt. trajectories)
u_lim = [CON.fmax; CON.phi_dot_max];