function x1 = oracle_MPC(x0, s, STR)

%--------------------------------------------------------------------------
% Context generation : 
% Using the horizon, a trajectory converging to the desired nominal
% trajectory is calculated.
% 
% Inputs : x0 - Initial Condition
%          s -  desired trajectory over the horizon 
%          STR - structure containing useful info
%
% Outputs : x1 - Next desired state
%
% Author: Okan Koc, IDSC ETH Zurich
% -------------------------------------------------------------------------

% release struct
bnds = STR.CGP.bounds;
dim_x = STR.dim_x;
horizon = size(s,2);

% set the bounds
x_bounds_left = -Inf(dim_x,horizon);
x_bounds_right = -x_bounds_left;
u_bounds_left = repmat(bnds(:,1),1,horizon);
u_bounds_right = repmat(bnds(:,2),1,horizon);
bounds_left = [x_bounds_left; u_bounds_left];
bounds_right = [x_bounds_right; u_bounds_right];
bounds(:,1) = bounds_left(:);
bounds(:,2) = bounds_right(:);

% CALL FMINCON
%%{
% fmincon has the advantage of starting from a nonzero value
eps = 0.1;
CxU = [s; u_bounds_left + eps * rand(size(u_bounds_left))];
cxu0 = CxU(:);
f = @(cxu) cost(cxu, x0, s, STR);
opts = optimset('Display', 'off', 'Algorithm','interior-point');
%opts = optimset('Display', 'off', 'Algorithm','sqp');
xmin = fmincon(f,cxu0,[],[],[],[],bounds(:,1),bounds(:,2),[],opts);
% only take the first step on the converging trajectory
x1 = xmin(1:dim_x);
%disp('FMINCON found xmin at:'); xmin
%}
end

function val = cost(testpt, x0, s, STR)

lambda = 1; % langrange-like multiplier

dim_x = STR.dim_x;
dim = STR.dim;
dims_tr = STR.dims_tr;
horizon = length(testpt)/(dim_x + dim);
% get the x's and u's out
cxu = reshape(testpt, dim_x + dim, horizon);
xs = cxu(1:dim_x,:);
xs = [x0, xs];
us = cxu(dim_x+1:end,:);

% add nominal trajectory error prediction to gp
% get necessary variables out
CON = STR.CON;
PAR = STR.PAR;
h = STR.h;
Q = STR.Q;
handle = STR.handle; % nominal function handle
tot_err = 0;

for i = 1:horizon
    pred = step_RK4(h,xs(:,i),us(:,i),CON,PAR,handle);
    err = pred - xs(:,i+1);
    tot_err = tot_err + lambda * err'*Q*err; % mu contribution
    dev = xs(dims_tr,i+1) - s(dims_tr,i);
    dev_err = dev' * dev; % trajectory deviation contribution
    tot_err = tot_err + dev_err;
end
val = tot_err;
% penalize input
%a = 1.5 * 1e-6;
%R = a * diag(length(testpt));
%pen = testpt'*R*testpt;
%val = nom_err + pen;

end