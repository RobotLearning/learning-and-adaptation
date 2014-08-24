function out = optim_trj(x_nom,x0,STR)

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
% form the difference matrix 
D = diag(ones(N-1,1),0) + diag(-ones(N-2,1),1);
D = [D, [zeros(N-2,1); -1]];
% vector for equating initial points
v = [1,zeros(1,N-1)];
% form the lipschitz vector punishing the derivatives of the opt. path
if nargin < 3
L = (x_nom(end) - x0) * ones(N-1,1); % straight line
else
L = L * ones(N-1,1); % TODO: does not have to be a constant vector
end

% call QP optimization with MATLAB
%%{
A = [D; -D];
b = [L; L];
% if older version remove '-convex'
opts = optimset('Display','off', 'Algorithm', 'interior-point-convex');
%opts = optimset('Display','iter', 'Algorithm', 'interior-point-convex');
x = quadprog(Q,-Q*x_nom,A,b,v,x0,[],[],[],opts);
%}

% call QP optimization with CVX
%{
cvx_begin
    % number of variables to optimize is N 
    l = 2; % 2,1,or Inf norms can be considered
    variable x(N) % converging trajectory
    minimize(norm(Q*(x-x_nom),l))
    subject to 
        D*x <= L %#ok
        D*x >= -L %#ok
        v*x == x0 %#ok
cvx_end
%}

% return the results
out = x;
% debug by plotting the results
%plot(1:N, x_nom, '-r', 1:N, x, '-b');