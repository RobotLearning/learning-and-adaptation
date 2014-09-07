function out = oracle_H(x_nom,x0,STR)

%--------------------------------------------------------------------------
% Context generation : 
% Using the horizon, a trajectory converging to the desired nominal
% trajectory is calculated.
% 
% Inputs : x_nom : desired next point on the trajectory
%          x0:     initial point, current state
%          STR:    containing bounds on the input u at time stage t
%                  contains also the Q matrix
%
% Author: Okan Koc, IDSC ETH Zurich
% -------------------------------------------------------------------------

% QP OPTIMIZATION
% optimize trajectory from x0 to endpoint 
N = length(x_nom);
% form the Q matrix
Q = STR.Q;
h = STR.h;

% TODO: replace for generic dynamics f 
bounds = STR.CGP.bounds;
Aeq = [1 0 0 0 0;
       0 0 1 0 0];
beq = x0([1,3]) + x0([2,4])*h;
A = [0 1 0 0 0;
     0 0 0 1 0;
     0 0 0 0 1];
A = [A; -A]; 
b = [max(x0(2) + h*bounds(1,1:2)'*sin(x0(5)));
     max(x0(4) + h*bounds(1,1:2)'*cos(x0(5)));
     max(x0(5) + h*bounds(2,1:2));
     -min(x0(2) + h*bounds(1,1:2)'*sin(x0(5)));
     -min(x0(4) + h*bounds(1,1:2)'*cos(x0(5)));
     -min(x0(5) + h*bounds(2,1:2));];
 
% adjust derivatives of x_nom
x_nom(2) = (x_nom(1) - x0(1))/h;
x_nom(4) = (x_nom(3) - x0(3))/h;


% if older version remove '-convex'
opts = optimset('Display','off', 'Algorithm', 'interior-point-convex');
%opts = optimset('Display','iter', 'Algorithm', 'interior-point-convex');
x = quadprog(Q,-Q*x_nom,A,b,Aeq,beq,[],[],[],opts);

% return the results
out = x;
% debug by plotting the results
%plot(1:N, x_nom, '-r', 1:N, x, '-b');