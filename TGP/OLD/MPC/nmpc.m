function xmin = nmpc(STR)

%--------------------------------------------------------------------------
% Nonlinear Model Predictive Control 
% Using DIRECT or FMINCON as optimization
%
% INPUTS:
% STR - structure containing useful info
%
% OUTPUTS:
% xmin - control input to be tried at the current iteration
%
% AUTHOR: Okan Koc, ETHZ IDSC Lab, 2012
%--------------------------------------------------------------------------

% release struct
if STR.FLAGS.BND
    bounds = STR.MPC.bounds;
else
    dim = STR.dim;
    bounds = Inf * [-ones(dim,1), ones(dim,1)];
end

% fmincon has the advantage of starting from a nonzero value
horizon = STR.MPC.HORIZON;
x0 = repmat(STR.MPC.u_past,horizon,1);
bounds = repmat(bounds,horizon,1);

% CALL DIRECT
%{
% Send options to Direct 
options.showits   = 0;
options.tol       = 0.05;
options.maxevals  = 200;
options.maxits    = 100;
options.maxdeep   = 100;
% Pass function as part of a Matlab Structure
problem.f = @(x) mpc(x, STR);
[~,xmin,hist] = Direct(problem,bounds,options); %#ok
%disp('DIRECT found xmin at:'); xmin
%figure(2); plot(hist(:,2),hist(:,3),'-*');
%}

% CALL FMINCON
%%{
f = @(x) mpc(x, STR);
%opts = optimset('Display', 'off', 'Algorithm','interior-point');
opts = optimset('Display', 'off', 'Algorithm','sqp');
xmin = fmincon(f,x0,[],[],[],[],bounds(:,1),bounds(:,2),[],opts);
% only take the first input at time t
xmin = xmin(1:length(STR.MPC.u_past));
%disp('FMINCON found xmin at:'); xmin
%}
end

function val = mpc(testpt, STR)

% add nominal trajectory error prediction to gp
% get necessary variables out
CON = STR.CON;
PAR = STR.PAR;
x_now = STR.MPC.x_now;
h = STR.h;
dim = STR.dim;
Q = STR.Q;
handle = STR.handle; % nominal function handle
horizon = STR.MPC.HORIZON;
trj = STR.MPC.TRJ;

% minimize the combined cost over the horizon
pred = zeros(size(trj));
nom_err = 0;
for i = 1:horizon
    pred(:,i) = step_RK4(h,x_now,testpt((dim*(i-1)+1):(dim*i)),CON,PAR,handle);
    x_now = pred(:,i);
    err = pred(:,i) - trj(:,i);
    nom_err = nom_err + err'*Q*err;
end
val = nom_err;
% penalize input
%a = 1.5 * 1e-6;
%R = a * diag(length(testpt));
%pen = testpt'*R*testpt;
%val = nom_err + pen;

end