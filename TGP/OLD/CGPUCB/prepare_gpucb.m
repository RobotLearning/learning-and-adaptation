%% PREPARE VARIABLES AND PARAMETERS FOR GPUCB
%--------------------------------------------------------------------------
% Called once by the GP main script (GP_Main)
%
% Author: Okan Koc, IDSC Lab, ETH Zurich
%
% -------------------------------------------------------------------------

% increase constraints for the learning process
if exist('increase_constr.m','file'), increase_constr; end

%%%%%%%%%%%%%%%%%%%%%%% GP-UCB CONSTANTS, PARAM, ETC. %%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% FUNCTION DRAWN FROM GP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct beta_t using Thm 2 in A.Krause'per
% hoeffding prob. constant
delta = 0.1;
% these constants depend on the kernel bound (see Krause's GP-UCB paper)
a = 1; b = 1;
% dimension of the cube
d = dim_ctx + dim;
% beta multiplier - more aggression necessary for tackling real problems
beta_mult = 0.2;
beta_f = @(t,delta,d,a,b,r) beta_mult * (2*log((t^2)*2*(pi^2)/(3*delta)) +...
       2*d*log((t^2)*d*b*r*sqrt(log(4*d*a/delta))));

% ALTERNATIVELY USE RKHS
% construct beta_f using Thm 3 in A.Krause's paper

%iteration/run numbers
iter = N-1;
% load into STR structures
STR.CON = CON;
STR.PAR = PAR;
STR.h = h;
STR.Q = Q;
STR.handle = fun_nom;
STR.dim = dim;
STR.dim_x = dim_x;
STR.dims_tr = dims_tr;
STR.CGP.HORIZON = horizon;
STR.FLAGS.BND = flag_bnd;
STR.FLAGS.DRV_BND = flag_drv_bnd;

% each run learns from the previous runs
runs = 3; 
err_total = zeros(runs);
% initialize past control signals
STR.CGP.u_past = zeros(dim,1);
STR.GPMPC.u_past = zeros(dim,1);
STR.MPC.u_past = zeros(dim,1);