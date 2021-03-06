function [u,STR] = cgp_ucb(t, cost, ctx, us, STR)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a function that generates a control signal to be tried
% INPUTS
% t - time stage of the last iteration
% costs - costs of the controls tried so far
% ctx - past contexts 
% us  - past controls tried
% STR - structure passing useful information
% OUTPUT
% u - control signal at time stage t
% STR - structure containing updated information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dim = STR.dim_act;
if t == 0, u = rand(dim,1); return; end

% past contexts + controls
cxu = [ctx; us];

meshsize = STR.meshsize;
var_noise = STR.var;
kern = STR.kernel;
lb = STR.lb;
ub = STR.ub;
Kinv = STR.Kinv; % previous inverse

%%%%%%%%%%%%%%%%% CONSTRUCT BETA CONSTANT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hoeffding prob. constant
delta = 0.1;
% bounds of the cube D
r = ub - lb;  
% these constants depends on the kernel bound (see Krause's GP-UCB paper)
a = 1; b = 1;
beta_mult = 0.05;
% Thm1 for discrete domain
beta_thm1 = @(t,D,delta) 2*log((t^2)*D*(pi^2)/(6*delta));
% Thm2 for continous compact domain
beta_thm2 = @(t,delta,d,a,b,r) beta_mult * (2*log((t^2)*2*(pi^2)/(3*delta)) +...
       2*d*log((t^2)*d*b*r*sqrt(log(4*d*a/delta))));
%beta = beta_thm2(t,delta,dim,a,b,r);
beta = beta_thm1(t,meshsize,delta);
STR.beta = beta;

%%%%%%%%%%%%%%%%% CONSTRUCT K MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform update for t
%%{
if t > 0
mat = ker_matrix_iter(t,cxu,kern,STR.K);
% add noise to last point
mat(end,end) = mat(end,end) + var_noise;
STR.K = mat;
end
%}

%%%%%%%%%%%%%%%%% CONSTRUCT INVERSE MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% perform update for t
%{
% calculate the matrix inverse incrementally
if (t > 1) % if t is 1, then Kinv is a scalar
% calculate A and store, calculate s
Q = ker_vector(t-1,cxu(t),cxu(1:t-1),kern);
A = Kinv * Q;
S = kern(cxu(t),cxu(t)) + var_noise; % usually 1
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
%%{
options = optimset('Display','off','Algorithm','interior-point');
options = optimset(options,'GradObj','off');
% fmincon has the advantage of starting from a nonzero value
u0 = rand(dim,1);
f = @(x) cgp(t,x,ctx,us,cost,STR);
u = fmincon(f,u0,[],[],[],[],lb,ub,[],options);
%disp('FMINCON found xmin at:'); u
%}

% CALL DIRECT
%{
% Send options to Direct 
options.showits   = 0;
options.tol       = 0.05;
options.maxevals  = 50;
options.maxits    = 50;
options.maxdeep   = 100;
% bounds
bounds(:,1) = lb;
bounds(:,2) = ub;
% Pass function as part of a Matlab Structure
problem.f = @(x) cgp(t,x,ctx,us,cost,STR);
[~,u,hist] = Direct(problem,bounds,options); %#ok
%disp('DIRECT found xmin at:'); u
%figure(2); plot(hist(:,2),hist(:,3),'-*');
%}

% CALL FMINCON WITH QP MINUS NORM
%%{
Q = STR.boost.Q;
[M,c_gp] = qp(t,ctx,us,cost,STR);
c = STR.boost.c + c_gp;
options = optimset('Display','off','Algorithm','sqp');
options = optimset(options,'GradObj','off');
% fmincon has the advantage of starting from a nonzero value
u0 = 0.5*ones(dim,1);
f = @(u) u'*Q*u + c'*u - sqrt(u'*M*u);
u = fmincon(f,u0,[],[],[],[],lb,ub,[],options);
%}

% QP IF 1D
%{
% form Q and c for QP 
Q_boost = STR.boost.Q;
c_boost = STR.boost.c;
Q_gp = 0; c_gp = 0;
if t > 0
[Q_gp, c_gp] = qp(t,ctx,us,cost,STR);
end
Q = Q_boost + Q_gp;
c = c_boost + c_gp;

%candidate for minimum - works in 1D case
u = -1/2 * (Q \ c);
if u > ub || u < lb
    val_ub = ub*Q*ub + ub*c;
    val_lb = lb*Q*lb + lb*c;
    if val_lb < val_ub
        u = lb;
    else
        u = ub;
    end
end
% ALTERNATIVE QP
% optimization toolbox should be installed
%options = optimset('Display','off','Algorithm','interior-point-convex');
%u = quadprog(2*Q,c,[],[],[],[],lb,ub,rand,options);
%}

end

%%%%%%%%%%%%%%% CONSTRUCT Q/c FOR OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Q,c] = qp(t,ctx,us,cost,STR)

dim = STR.dim_act;
beta = STR.beta;
mat = STR.K;
%mat = STR.Kinv;
sigma_s2 = STR.sigma_s2;
ctx_now = STR.ctx;
ker_ctx = STR.ker_ctx;
kappa = ker_vector(t,ctx_now,ctx,ker_ctx);
Kappa = repmat(kappa,1,dim) .* us';
L = STR.l_act;
InvLambda = diag(1./(L.^2));

c = InvLambda * Kappa' * (mat \ cost(:));
if dim == 1
    normu2 = beta * (sigma_s2 * InvLambda - ...
                InvLambda * Kappa' * (mat \ Kappa) * InvLambda);
    c_sigma = sqrt(normu2);              
    c = c - c_sigma;
    Q = 0;
else
    Q = beta * (sigma_s2 * InvLambda - ...
                InvLambda * Kappa' * (mat \ Kappa) * InvLambda);
end

end

%%%%%%%%%%%% BOOSTED FUNCTION TO BE OPTIMIZED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = cgp(t, u, ctx, us, cost, STR)

ctx_now = STR.ctx; % context now
cxu = [ctx; us];
boost = STR.boost.fun;
pt = [ctx_now; u];
boost_val = boost(u);
% call gp after getting one cost (for conditioning)
if t > 0
    ucb_val = mu_beta_sigma(t,pt,cxu,cost,STR);
    val = boost_val + ucb_val;
else
    val = boost_val;
end

% TAKING MAXIMUM OF A BUNCH OF TEST POINTS
%{
testsize = 100;
xx = linspace(lb, ub, testsize);
xxc = [repmat(ctx_now, 1, testsize); xx];
[m,s2,~] = mu_and_sigma_vec(t,xxc,cxu,cost,STR);
vals = m - sqrt(beta * s2);
[~,ind] = min(vals);
u = xx(ind);
%}
end

%%%%%%%%%%%% CONSTRUCT UCB FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [val,der] = mu_beta_sigma(t,x,cxu,cost,STR)

% release STR
dim = STR.dim;
kern = STR.kernel;
mat = STR.K;
%invmat = STR.Kinv;
beta = STR.beta;
l = STR.l;
sigma_s2 = STR.sigma_s2;
% calculate kernel vector k and store
ker_vec = ker_vector(t,x,cxu,kern);
% create mu function 
mu = @(x) ker_vec' * (mat \ cost(:));
%mu = @(x) ker_vec' * invmat * cost(:);

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
rep = (xs - cxu) .* ker_vec_rep'; 
der = -sigma_s2*(diag(l)^2) \ (rep * (mat \ y));
%der = -sigma_s2*(diag(l)^2) \ (rep * invmat * y);

end

%%%%%%%%%%%% CONSTRUCT UCB VECTOR-VALUED FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: does it work for cgp-ucb?
% here x is actually a matrix of inputs
function [m,s2,ders] = mu_and_sigma_vec(t,x,cxu,cost,STR)

% release STR
kern = STR.kernel;
mat = STR.K;
%invmat = STR.Kinv;
l = STR.l;
% calculate kernel vector k and store
ker_mat = zeros(size(cxu,2), size(x,2));
for i = 1:size(x,2)
    ker_mat(:,i) = ker_vector(t,x(:,i),cxu,kern);
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
