% GPUCB optimization using sampling
% 1-dimensional, no context!

function [ind, u, Kinv] = gpucb(t, cost, u_past, Kinv, kern)

% if at start point return one of the endpoints
 if t == 0, ind = 1; u = 0; return; end

% optimization assumptions
var_noise = 0.025;
lb = 0;
ub = 1;
 
% CONSTRUCT BETA CONSTANT 
% hoeffding prob. constant
delta = 0.1;
% size of the mesh
meshsize = 100;
% Thm1 for discrete space
beta_thm1 = @(t,D,delta) 2*log((t^2)*D*(pi^2)/(6*delta));
beta = beta_thm1(t,meshsize,delta);

% CONSTRUCT INVERSE MATRIX 
% perform update for t
% calculate the matrix inverse incrementally
if (t > 1) % if t is 1, then Kinv is a scalar
    % calculate A and store, calculate s
    Q = ker_vector(t-1,u_past(t),u_past(1:t-1),kern);
    A = Kinv * Q;
    S = kern(u_past(t), u_past(t)) + var_noise; % usually 1
    % calculate M
    m = 1/(S - Q'*A);
    % calculate Phat
    P_hat = Kinv + m*(A*A');
    % calculate Qhat
    Q_hat = -m*A;
    % form the full inverse matrix
    Kinv = [P_hat, Q_hat; Q_hat', m];
end

% TAKING MAXIMUM OF A BUNCH OF TEST POINTS
xx = linspace(lb, ub, meshsize);
[m,s2] = mu_and_sigma_vec(t,xx,u_past,cost,Kinv,kern);
vals = m + sqrt(beta * s2);
[~,ind] = max(vals);
u = xx(ind);

end

%%%%%%%%%%%% CONSTRUCT UCB VECTOR-VALUED FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%

% here x is actually a matrix of inputs
function [m,s2] = mu_and_sigma_vec(t,x,u_past,cost,Kinv,kern)

%mat = K;
invmat = Kinv;
% calculate kernel vector k and store
ker_mat = zeros(size(u_past,2), size(x,2));
for i = 1:size(x,2)
    ker_mat(:,i) = ker_vector(t,x(:,i),u_past,kern);
end    
% calculate prior covariance
ker_prior_mat = ker_matrix(x,kern);
% create mu function 
%mu = @(x) ker_mat' * (mat \ cost(:));
mu = @(x) ker_mat' * invmat * cost(:);

% create sigma function
% this makes sure we don't have negative values due to small num. errors
%cov = @(x) max(ker_prior_mat - ker_mat' * (mat \ ker_mat),0);
cov = @(x) max(ker_prior_mat - ker_mat' * invmat * ker_mat, 0);

% pass mu and sigma functions as function handles
m = mu(x);
s2 = diag(cov(x));

end
