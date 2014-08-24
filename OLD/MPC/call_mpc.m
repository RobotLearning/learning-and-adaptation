% used for the horizon in NMPC
%STR.MPC.TRJ = ctx(:,2:end);
STR.MPC.TRJ = x_hor(:,t+1:t+horizon);
STR.MPC.x_now = x_MPC(:,t);

% call the mpc code
u_mpc(:,t) = nmpc(STR);

% simulate one-step with a competitive method (NMPC)
x_MPC(:,t+1) = step_RK4(h,x_MPC(:,t),u_mpc(:,t),CON,PAR,fun_real);

% for fmincon to start from a nonzero value
STR.MPC.u_past = u_mpc(:,t); 

% update bounds [using derivative constraints]
if STR.FLAGS.DRV_BND
STR.MPC.bounds = update_bnds(t,STR,u_mpc(:,t));
end