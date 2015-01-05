% Sets system constants, constraints, disturbance model

%System constants
par.Iy = 0.0023;
par.L = 0.17;
par.m = 0.468;
par.g = 9.81;            
par.Quad.A = 0.4 * 0.4; % 40 cm by 40 cm quadrocopter
            
%Single thrust constraints
con.fi_max = 4.5;
con.fi_min = 0.4;
con.fi_dot_max = 27;
%Collective thrust constraints
con.fmax = 4 * con.fi_max;
con.fmin = 4 * con.fi_min;
con.f_dot_max = 4 * con.fi_dot_max;
%Angle constraints
con.phi_max = pi/2; % only used for trajectory generation
con.phi_dot_max = 22;
con.phi_ddot_max = 150;
%Number of constraint evaluation points
con.NumT = 500;

%cost matrix
% Scaling matrix Sw
Sw = diag([1, 1, 1, 1, 0.1]); 
% Q matrix appears in the cost function
Q = Sw'*Sw;

dist.Type = 'gravity';
dist.g = 10.5; % in Uranus
%dist.g = 3.711; % in Mars

%{
dist.Type = 'wind'; % disturbance type (or actuator)
dist.Angle = 0; % angle of the blowing disturbance
dist.BlowPressure = 60; % Fan blow pressure N/m2
%}