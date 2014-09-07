%% dynamics_osci.m
% *Summary:* Implements ths ODE for simulating the oscilator dynamics, where 
% an input torque f can be applied  
%
%    function dz = dynamics_osci(t,z,u)
%    
%    see Nonlinear Optimal Control: A Control Lyapunov Function and 
%       Receding Horizon Perspective for details on the derivation
%
% *Input arguments:*
%
%	t     current time step (called from ODE solver)
%   z     state                                                    [2 x 1]
%   u     (optional): torque f(t) applied to pendulum
%
% *Output arguments:*
%   
%   dz    if 3 input arguments:      state derivative wrt time
%
%   Note: It is assumed that the state variables are of the following order:
%         dtheta:  [rad/s] angular velocity of pendulum
%         theta:   [rad]   angle of pendulum
%

function dz = dynamics_osci(t,z,u)
% function dz = dynamics_osci(t,z) % this is for testing the dynamics

% u = 3*z(2); % this is the optimal input w.r.t. \int_0^{+\infty}z^2(2)+u^2

dz = zeros(2,1);
dz(1) = z(2);
dz(2) = -z(1)*(pi/2 + atan(5*z(1))) ...
        - 2.5*z(1)*z(1)/(1+25*z(1)*z(1)) ...
        + 4*z(2) ...
        + 3*u;