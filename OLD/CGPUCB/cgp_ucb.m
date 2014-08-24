%% CGP-UCB code
%--------------------------------------------------------------------------
% Called at each iteration once by the GP main script (GP_Main)
%
% Author: Okan Koc, IDSC Lab, ETH Zurich
%
% -------------------------------------------------------------------------

%x1 = oracle(x_hor(:,t+1),x_GP(:,t),STR);
x1 = x_hor(:,t+1);
STR.ctx = [x_GP(:,t); x1];

% get boost matrix Q and vector c
ufake = zeros(dim,1);
[~, aout, Bout] = fun_nom(t,x_GP(:,t),PAR,ufake,true);
STR.boost.Q = (h^2)*Bout'*Q*Bout;
STR.boost.c = 2*h*Bout'*Q*(x_GP(:,t) + h*aout - x1);

% number of data points GP is conditioned on
num_data = (run-1) * iter + t - 1;

% beta is slightly increasing at each iteration
%STR.beta = beta_f(t,delta,d,a,b,r);
STR.beta = beta_f(num_data + 1,delta,d,a,b,r);

% find new control signal to be tried
% conditioning on the previous data is determined by past var.
% condition on previous runs
past = past_example_data + (1: window);
window = window + 1;
STR.CGP.x_now = x_GP(:,t);

if flag_GPMPC
    STR.GPMPC.TRJ = x_hor(:,t+1:t+horizon);
    u(:,t) = gp_mpc(GPSTR, CxU(past,:), cost(past), STR);
    STR.GPMPC.u_past = u(:,t);
else
    u(:,t) = gp_ucb_boost_R(GPSTR, CxU(past,:),cost(past),STR);
    phi = STR.ctx(5);
    dy_dot = STR.ctx(2) - STR.ctx(7);
    dbg(t) = STR.DEBUG.FNC(u(:,t),phi,dy_dot);
    %{
    ctxs = CxU(past,1:dim_ctx)';
    us = CxU(past,dim_ctx+1:end)';
    [u(:,t),STR] = gp_ucb_boost(num_data,cost(past),ctxs,us,STR,GPSTR);
    %}
    
    % for fmincon to start from a nonzero value
    STR.CGP.u_past = u(:,t);
end

% update bounds [using derivative constraints]
if STR.FLAGS.DRV_BND
STR.CGP.bounds = update_bnds(t,STR,u(:,t));
end

% simulate one-step with GP
x_GP(:,t+1) = step_RK4(h,x_GP(:,t),u(:,t),CON,PAR,fun_real);

%%%%%%%%%%%%%%%%% ADD NOISE %%%%%%%%%%%%%%%%%%%
%x_GP(:,t+1) = x_GP(:,t+1) + sqrt(eps) * randn;

dev = x_GP(:,t+1) - x1;
newCost = dev'*Q*dev;
% subtract nominal model prediction error
pred = step_RK4(h,x_GP(:,t),u(:,t),CON,PAR,fun_nom);
pred_dev = pred - x1;
pred_err = pred_dev'*Q*pred_dev;
newCost = newCost - pred_err;
% DEBUG!!!!
%newCost = newCost - 2*h*8*(x_GP(2,t)-x_hor(2,t+1))*sin(x_GP(5,t));

% add u and newCost to training data/observations pair
present = past_example_data + window;
CxU(present,:) = [STR.ctx; u(:,t)]';
cost(present) = newCost;