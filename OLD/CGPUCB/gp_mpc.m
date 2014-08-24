function xmin = gp_mpc(gpstr, cxu, y, STR)

%--------------------------------------------------------------------------
% CGP-MPC 
% roll-out algorithm using learned dynamics
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
%-----------------

% release struct
hyp = gpstr.hyp;
inf = gpstr.inf;
meanfunc = gpstr.meanfunc;
covfunc = gpstr.covfunc;
likfunc = gpstr.likfunc;

% are we bounding the search space?
if STR.FLAGS.BND
    bounds = STR.CGP.bounds;
else
    dim = STR.dim;
    bounds = Inf * [-ones(dim,1), ones(dim,1)];
end

horizon = STR.GPMPC.HORIZON;
x0 = repmat(STR.GPMPC.u_past,horizon,1);
bounds = repmat(bounds,horizon,1);

f = @(x) gpmpc(hyp, inf, meanfunc, covfunc, likfunc, cxu, y, x, STR);
opts = optimset('Display', 'off', 'Algorithm','interior-point');
xmin = fmincon(f,x0,[],[],[],[],bounds(:,1),bounds(:,2),[],opts);
xmin = xmin(1:length(STR.GPMPC.u_past));

end

function val = gpmpc(hyp, inf, meanfunc, covfunc, likfunc, ...
                    u, y, testpt, STR)

% add nominal trajectory error prediction to gp
% get necessary variables out
CON = STR.CON;
PAR = STR.PAR;
x_now = STR.CGP.x_now;
h = STR.h;
dim = STR.dim;
Q = STR.Q;
handle = STR.handle; % nominal function handle
horizon = STR.GPMPC.HORIZON;
trj = STR.GPMPC.TRJ;

% minimize the combined cost over the horizon
pred = zeros(size(trj));
nom_err = 0;
for i = 1:horizon
    ctx = [x_now; trj(:,i)];
    pt = [ctx; testpt((dim*(i-1)+1):(dim*i))]';
    pred(:,i) = step_RK4(h,x_now,testpt((dim*(i-1)+1):(dim*i)),CON,PAR,handle);
    x_now = pred(:,i);
    err = pred(:,i) - trj(:,i);
    [m, ~] = gp(hyp, inf, meanfunc, covfunc, likfunc, u, y, pt);
    err = err'*Q*err + m;
    nom_err = nom_err + err;
end
val = nom_err;

end