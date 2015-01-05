function [u,STR] = gp_ucb_boost(t, cost, ctx, us, STR, GPSTR)

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

%%%%%%%%%%%%%%%%%% FORMAT DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim_ctx = STR.dim_x * 2;
dim = STR.dim;
%if t == 0, u = rand(dim,1); return; end
% extract bounds for the optimization
if STR.FLAGS.BND
    bounds = STR.CGP.bounds;
else
    bounds = Inf * [-ones(dim,1), ones(dim,1)];
end

hyp_ctx = Inf(1,dim_ctx);
maskCtx = GPSTR.covfunc{2}{1}{2}{1};
ctx_l = sum(maskCtx);
hyp_ctx(maskCtx > 0) = GPSTR.hyp.cov(1:ctx_l);
% if ard
%hyp_act = Inf(1,dim);
%maskAct = GPSTR.covfunc{2}{2}{2}{1};
%hyp_act(find(maskAct(dim_ctx+1:end))) = GPSTR.hyp.cov(ctx_l+2:end);
hyp_act = [0,Inf];

%%%%%%%%%%%%%%%%%% RELEASE HYPERPARAMETERS AND CONSTRUCT KERNEL %%%%%%%%%%%
% variance of the estimated noise
var_noise = exp(GPSTR.hyp.lik)^2;
% kernel function for the contexts
ctx_hyp = exp(hyp_ctx);
ctx_scale = exp(GPSTR.hyp.cov(ctx_l+1))^2;
ctx_type = 'squared exponential ard';
ctx_kern = @(x1,x2) ctx_scale * kernel(x1(1:dim_ctx), ...
                                       x2(1:dim_ctx), ctx_hyp, ctx_type);
% kernel function for the actions
act_hyp = []; %exp(hyp_act);
act_scale = 1; %exp(hyp.cov(end));
act_type = 'linear iso';
act_kern = @(x1,x2) act_scale * kernel(x1(dim_ctx+1), ...
                                       x2(dim_ctx+1), ...
                                       act_hyp, act_type);
kern = @(x1,x2) act_kern(x1,x2) * ctx_kern(x1,x2);
STR.L = act_hyp;
STR.ctx_scale = ctx_scale;
STR.ctx_kern = ctx_kern;
STR.kernel = kern;

% past contexts + controls
cxu = [ctx; us];
% previous inverse
Kinv = STR.Kinv; 

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
Q = ker_vector(t-1,cxu(:,t),cxu(:,1:t-1),kern);
A = Kinv * Q;
S = kern(cxu(:,t),cxu(:,t)) + var_noise; % usually 1
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

% debug m and s2 by plotting
%plotsurface(t, cxu, cost, STR);

%%%%%%%%%%%%%%%%%%% FUNCTION OPTIMIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CALL FMINCON
%%{
f = @(x) gpucb(t, cxu, cost, x, STR);
%opts = optimset('Display', 'off', 'Algorithm','interior-point');
opts = optimset('Display','off','Algorithm', 'sqp');
% fmincon has the advantage of starting from a nonzero value
u0 = STR.CGP.u_past;
u = fmincon(f,u0,[],[],[],[],bounds(:,1),bounds(:,2),[],opts);
%disp('FMINCON found xmin at:'); xmin
%}

% CALL FMINCON WITH QP MINUS NORM
%{
Q = STR.boost.Q;
[M,c_gp] = qp(t,ctx,us,cost,STR);
c = STR.boost.c + c_gp;
%options = optimset('Display','off','Algorithm','interior-point');
options = optimset('Display','off','Algorithm', 'sqp');
options = optimset(options,'GradObj','off');
% fmincon has the advantage of starting from a nonzero value
u0 = STR.CGP.u_past;
f = @(u) u'*Q*u + c'*u - sqrt(u'*M*u);
u = fmincon(f,u0,[],[],[],[],bounds(:,1),bounds(:,2),[],options);
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

dim = STR.dim;
beta = STR.beta;
mat = STR.K;
%mat = STR.Kinv;
ctx_scale = STR.ctx_scale;
ctx_now = STR.ctx;
ker_ctx = STR.ctx_kern;
kappa = ker_vector(t,ctx_now,ctx,ker_ctx);
Kappa = repmat(kappa,1,dim) .* us';
L = STR.L;
InvLambda = diag(1./(L.^2));

c = InvLambda * Kappa' * (mat \ cost(:));
%c = InvLambda * Kappa' * mat * cost(:);
if dim == 1
    normu2 = beta * (ctx_scale * InvLambda - ...
                InvLambda * Kappa' * (mat \ Kappa) * InvLambda);
    %normu2 = beta * (ctx_scale * InvLambda - ...
    %            InvLambda * Kappa' * mat * Kappa * InvLambda);
    c_sigma = sqrt(normu2);              
    c = c - c_sigma;
    Q = 0;
else
    Q = beta * (ctx_scale * InvLambda - ...
                InvLambda * Kappa' * (mat \ Kappa) * InvLambda);
    %Q = beta * (ctx_scale * InvLambda - ...
    %            InvLambda * Kappa' * mat * Kappa * InvLambda);
end

end

%%%%%%%%%%%% BOOSTED FUNCTION TO BE OPTIMIZED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = gpucb(t, cxu, cost, u, STR)

% add nominal trajectory error prediction to gp
% get necessary variables out
CON = STR.CON;
PAR = STR.PAR;
x_now = STR.CGP.x_now;
h = STR.h;
ctx = STR.ctx;
Q = STR.Q;
dim_x = length(x_now);
beta = STR.beta;
kern = STR.kernel;
%mat = STR.Kinv;
mat = STR.K;
handle = STR.handle; % nominal function handle
pred = step_RK4(h,x_now,u,CON,PAR,handle);

% minimize the combined cost
nom_err = pred(:) - ctx(dim_x + (1:dim_x));
nom_err = nom_err'*Q*nom_err;

% call gp after getting one cost (for conditioning)
if t > 0
pt = [ctx; u]';
ucb_val = mu_beta_sigma(t,pt,cxu,cost,beta,kern,mat);
val = nom_err + ucb_val;
else
val = nom_err;
end

% debug beta level (exploration)
%fprintf('nom = %f, m = %f, s2 = %f, explore = %f.\n', ...
%        nom_err, m, s2, sqrt(beta * s2));

end

%%%%%%%%%%%% UCB FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function val = mu_beta_sigma(t,x,cxu,cost,beta,kern,mat)

% calculate kernel vector k and store
vec = ker_vector(t,x,cxu,kern);
% create mu function 
mu = vec' * (mat \ cost(:));
%mu = vec' * mat * cost(:);
% create sigma function
s2 = max(kern(x,x) - vec' * (mat \ vec), 0);
%s2 = max(kern(x,x) - vec' * mat * vec, 0);
% pass mu and sigma functions as function handles
val = mu - sqrt(beta * s2);

end

%%%%%%%%%%%% CONSTRUCT UCB VECTOR-VALUED FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%

% here x is actually a matrix of inputs
function [m,s2] = mu_and_sigma_vec(t,x,u_past,cost,STR)

% release STR
kern = STR.kernel;
mat = STR.K;
%invmat = STR.Kinv;
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

end

function plotsurface(t, cxu, y, STR)

bounds = STR.CGP.bounds;
CON = STR.CON;
PAR = STR.PAR;
x_now = STR.CGP.x_now;
h = STR.h;
ctx = STR.ctx;
Q = STR.Q;
dim_x = length(x_now);
beta = STR.beta;
handle = STR.handle; % nominal function handle
mesh1 = bounds(1,1):1:bounds(1,2);
mesh2 = bounds(2,1):1:bounds(2,2);
meshx = repmat(mesh1, length(mesh2), 1);
meshy = repmat(mesh2, 1, length(mesh1));
meshxx = meshx(:);
meshyy = meshy(:);
X = [meshxx, meshyy];

% get the nominal boost (mean)
mboost = zeros(length(X),1);
for i = 1:length(X)
    pred = step_RK4(h,x_now,X(i,:)',CON,PAR,handle);
    % minimize the combined cost
    err = pred(:) - ctx(dim_x + (1:dim_x));
    mboost(i) = err'*Q*err;
end

ctx = STR.ctx;
CTX = repmat(ctx', length(X), 1);
PTS = [CTX, X];
[M,S2] = mu_and_sigma_vec(t,PTS',cxu,y,STR);
[U1s,U2s] = meshgrid(mesh1,mesh2);
Boost = reshape(mboost, length(mesh2), length(mesh1));
Ms = reshape(M, length(mesh2), length(mesh1));
Ss = reshape(S2, length(mesh2), length(mesh1));
figure(4); surf(U1s, U2s, Boost);
figure(5); surf(U1s, U2s, Ms);
figure(6); surf(U1s, U2s, sqrt(beta*Ss));

end
