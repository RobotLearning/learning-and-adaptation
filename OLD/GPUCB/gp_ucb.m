function [u,STR] = gp_ucb(t, cost, u_past, STR)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a function that generates a control signal to be tried
% INPUTS
% t - time stage of the last iteration
% costs - costs of the controls tried so far
% u_past - past control signals tried
% STR - structure passing useful information
% OUTPUT
% u - control signal at time stage t
% STR - structure containing updated information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if at start point return one of the endpoints
 if t == 0, u = 0; return; end

dim = STR.dim;
meshsize = STR.meshsize;
var_noise = STR.var;
kern = STR.kernel;
lb = STR.lb;
ub = STR.ub;
%Kinv = STR.Kinv; % previous inverse

%%%%%%%%%%%%%%%%% CONSTRUCT BETA CONSTANT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hoeffding prob. constant
delta = 0.1;
% bounds of the cube D
r = ub - lb;  
% these constants depends on the kernel bound (see Krause's GP-UCB paper)
a = 1; b = 1;
beta_mult = 0.05;
% Thm1 for discrete space
beta_thm1 = @(t,D,delta) 2*log((t^2)*D*(pi^2)/(6*delta));
% Thm2 for continous compact space
beta_thm2 = @(t,delta,d,a,b,r) beta_mult * (2*log((t^2)*2*(pi^2)/(3*delta)) +...
       2*d*log((t^2)*d*b*r*sqrt(log(4*d*a/delta))));
%beta = beta_thm2(t,delta,dim,a,b,r);
beta = beta_thm1(t,meshsize,delta);
STR.beta = beta;

%%%%%%%%%%%%%%%%% CONSTRUCT K MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform update for t
%%{
mat = ker_matrix_iter(t,u_past,kern,STR.K);
% add noise to last point
mat(end,end) = mat(end,end) + var_noise;
STR.K = mat;
%}

%%%%%%%%%%%%%%%%% CONSTRUCT INVERSE MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform update for t
%{
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
STR.Kinv = Kinv;
end
%}

%%%%%%%%%%%%%%%%%%% FUNCTION OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALL FMINCON
%{
options = optimset('Display','off','Algorithm','interior-point');
options = optimset(options,'GradObj','off');
% fmincon has the advantage of starting from a nonzero value
u0 = 1;
f = @(x) mu_beta_sigma(t,x,u_past,cost,STR);
u = fmincon(f,u0,[],[],[],[],lb,ub,[],options);
%disp('FMINCON found xmin at:'); u
%}

% CALL DIRECT
%%{
% Send options to Direct 
options.showits   = 0;
options.tol       = 0.05;
options.maxevals  = 100;
options.maxits    = 100;
options.maxdeep   = 100;
% bounds
bounds(:,1) = lb;
bounds(:,2) = ub;
% Pass function as part of a Matlab Structure
problem.f = @(x) mu_beta_sigma(t,x,u_past,cost,STR);
[~,u,hist] = Direct(problem,bounds,options); %#ok
%disp('DIRECT found xmin at:'); u
%figure(2); plot(hist(:,2),hist(:,3),'-*');
%}

% TAKING MAXIMUM OF A BUNCH OF TEST POINTS
%{
testsize = 100;
xx = linspace(lb, ub, testsize);
[m,s2,~] = mu_and_sigma_vec(t,xx,u_past,cost,STR);
vals = m - sqrt(beta * s2);
[~,ind] = min(vals);
u = xx(ind);
%}

% PLOTTING GP-UCB ACTION
%{
figure(3);
testsize = 100;
xx = linspace(lb, ub, testsize);
[m,s2,~] = mu_and_sigma_vec(t,xx,u_past,cost,STR);
f = [m+sqrt(beta)*sqrt(s2); flipdim(m-sqrt(beta)*sqrt(s2),1)]; 
fill([xx(:); flipdim(xx(:),1)], f, [7 7 7]/8)
hold on; 
plot(xx, m, '-r'); 
plot(u, mu_beta_sigma(t,u,u_past,cost,STR), '+');
hold off;
%}

end

%%%%%%%%%%%% CONSTRUCT UCB FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [val,der] = mu_beta_sigma(t,x,u_past,cost,STR)

% release STR
dim = STR.dim;
kern = STR.kernel;
mat = STR.K;
%invmat = STR.Kinv;
beta = STR.beta;
l = STR.l;
sigma_s2 = STR.sigma_s2;
% calculate kernel vector k and store
ker_vec = ker_vector(t,x,u_past,kern);
% create mu function 
mu = @(x) ker_vec' * (mat \ cost(:));
%mu = @(x) ker_mat' * invmat * cost(:);

% create sigma function
% this makes sure we don't have negative values due to small num. errors
sigma2 = @(x) max(kern(x,x) - ker_vec' * (mat \ ker_vec),0);
%sigma2 = @(x) max(kern(x,x) - ker_vec' * invmat * ker_vec, 0);

% pass mu and sigma functions as function handles
val = mu(x) - sqrt(beta * sigma2(x));

% derivative assuming exponential kernel!
y = cost(:) + 2 * sqrt(beta) * ker_vec;
xs = repmat(x,1,t);
ker_vec_rep = repmat(ker_vec,1,dim);
rep = (xs - u_past) .* ker_vec_rep'; 
der = -sigma_s2*(diag(l)^2) \ (rep * (mat \ y));

end

%%%%%%%%%%%% CONSTRUCT UCB VECTOR-VALUED FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%

% here x is actually a matrix of inputs
function [m,s2,ders] = mu_and_sigma_vec(t,x,u_past,cost,STR)

% release STR
kern = STR.kernel;
mat = STR.K;
%invmat = STR.Kinv;
l = STR.l;
% calculate kernel vector k and store
ker_mat = zeros(size(u_past,2), size(x,2));
for i = 1:size(x,2)
    ker_mat(:,i) = ker_vector(t,x(:,i),u_past,kern);
end    
% calculate prior covariance
ker_prior_mat = ker_matrix(x,kern);
% create mu function 
mu = @(x) ker_mat' * (mat \ cost(:));
%mu = @(x) ker_mat' * invmat * cost(:);

% create sigma function
% this makes sure we don't have negative values due to small num. errors
cov = @(x) max(ker_prior_mat - ker_mat' * (mat \ ker_mat),0);
%cov = @(x) max(ker_prior_mat - ker_mat' * invmat * ker_mat, 0);

% pass mu and sigma functions as function handles
m = mu(x);
s2 = diag(cov(x));

ders = 0;
% TODO: update ders for multi-D, multi-test point case
% derivative assuming exponential kernel!
%y = cost(:) + 2 * sqrt(beta) * ker_mat;
%ders = -1/(l^2) * ((x - u_past(:)) .* ker_mat)' * (mat \ y);

end
