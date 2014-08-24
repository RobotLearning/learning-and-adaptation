function xmin = gp_ucb_boost_R(gpstr, cxu, y, STR)

%--------------------------------------------------------------------------
% CGP-UCB boosted with nominal model (i.e. providing KNOWN mean function)
% Using DIRECT or FMINCON as optimization
% THIS VERSION (COMPARED TO GP_UCB_BOOST) USES RASMUSSEN'S GP FUNCTION
%
% INPUTS:
% cxu - past contexts + past control inputs tried
% y - costs incurred in the past time stages
% gpstr - (hyperparameter) structure learned using RASMUSSSEN's code
% STR - structure containing useful info
%
% OUTPUTS:
% xmin - control input to be tried at the current iteration
%
% AUTHOR: Okan Koc, ETHZ IDSC Lab, 2012
%--------------------------------------------------------------------------

% release struct
hyp = gpstr.hyp;
inf = gpstr.inf;
meanfunc = gpstr.meanfunc;
covfunc = gpstr.covfunc;
likfunc = gpstr.likfunc;
beta = STR.beta;

% are we bounding the search space?
if STR.FLAGS.BND
    bounds = STR.CGP.bounds;
else
    dim = STR.dim;
    bounds = Inf * [-ones(dim,1), ones(dim,1)];
end

% debug m and s2 by plotting
%plotsurface(hyp, inf, meanfunc, covfunc, likfunc, cxu, y, beta, STR);

% fmincon has the advantage of starting from a nonzero value
x0 = STR.CGP.u_past;

% CALL DIRECT
%{
% Send options to Direct 
options.showits   = 0;
options.tol       = 0.05;
options.maxevals  = 200;
options.maxits    = 100;
options.maxdeep   = 100;
% Pass function as part of a Matlab Structure
problem.f = @(x) gpucb(hyp, inf, meanfunc, covfunc, likfunc, cxu, y, x, beta, STR);
[~,xmin,hist] = Direct(problem,bounds,options); %#ok
%disp('DIRECT found xmin at:'); xmin
%figure(2); plot(hist(:,2),hist(:,3),'-*');
%}

% CALL FMINCON
%%{
f = @(x) gpucb(hyp, inf, meanfunc, covfunc, likfunc, cxu, y, x, beta, STR);
%opts = optimset('Display', 'off', 'Algorithm','interior-point');
opts = optimset('Display', 'off', 'Algorithm','sqp');
%disp('FMINCON found xmin at:'); xmin
xmin = fmincon(f,x0,[],[],[],[],bounds(:,1),bounds(:,2),[],opts);
%}

end

function val = gpucb(hyp, inf, meanfunc, covfunc, likfunc, ...
                    u, y, testpt, beta, STR)

% add nominal trajectory error prediction to gp
% get necessary variables out
CON = STR.CON;
PAR = STR.PAR;
x_now = STR.CGP.x_now;
h = STR.h;
ctx = STR.ctx;
Q = STR.Q;
dim_x = length(x_now);
handle = STR.handle; % nominal function handle

% get the nominal boost (mean)
pred = step_RK4(h,x_now,testpt,CON,PAR,handle);
% minimize the combined cost
nom_err = pred(:) - ctx(dim_x + (1:dim_x));
nom_err = nom_err'*Q*nom_err;
% call gp
pt = [ctx; testpt]';
[m, s2] = gp(hyp, inf, meanfunc, covfunc, likfunc, u, y, pt);
val = nom_err + m - sqrt(beta * s2);

%{
% debug function
% assuming weight matrix is eye(6)
phi = ctx(5);
dy_dot = ctx(2) - ctx(7);
u = testpt;
dbg_cost = debug(u,phi,dy_dot);
val = nom_err + dbg_cost;
%}

% debug beta level (exploration)
%fprintf('nom = %f, m = %f, s2 = %f, explore = %f.\n', ...
%        nom_err, m, s2, sqrt(beta * s2));

end

function plotsurface(hyp, inf, meanfunc, covfunc, likfunc, ...
                        cxu, y, beta, STR)

bounds = STR.CGP.bounds;
CON = STR.CON;
PAR = STR.PAR;
x_now = STR.CGP.x_now;
h = STR.h;
ctx = STR.ctx;
Q = STR.Q;
dim_x = length(x_now);
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
[M,S2] = gp(hyp,inf,meanfunc,covfunc,likfunc,cxu,y,PTS);
[U1s,U2s] = meshgrid(mesh1,mesh2);
Boost = reshape(mboost, length(mesh2), length(mesh1));
Ms = reshape(M, length(mesh2), length(mesh1));
Ss = reshape(S2, length(mesh2), length(mesh1));
close all;
figure(1); surf(U1s, U2s, Boost);
figure(2); surf(U1s, U2s, Ms);
figure(3); surf(U1s, U2s, sqrt(beta*Ss));

end