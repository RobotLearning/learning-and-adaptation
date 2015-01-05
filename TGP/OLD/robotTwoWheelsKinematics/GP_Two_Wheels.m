%% Set system parameters

clc; clear; close all;
% run ILC n many iterations and get n many errors
% u_app and err are defined in the ILC-script
if ~exist('GP_Two_Wheels.mat','file')
    ILC_Two_Wheels;
    save('GP_Two_Wheels.mat','t', 'x', 'err','u_app','STR');
    clc; clear; close all;
end
load('GP_Two_Wheels.mat'); 
% extract necessary variables from structure
CON = STR.CON;
PAR = STR.PAR;
COV = STR.COV;
% remove the other variables (used in ILC) to avoid overhead
STR = rmfield(STR, ... 
{'F','G','H','C','umin','umax','x_con','d0','P0','S','w','lsq_cost'});
% u is constrained to be less than max. allowed
u_cnstr = min(CON.w1_cnstr,CON.w2_cnstr);
% number of iterations
% n = size(err,2);
% only use 10 of the examples
n = 10;
% dimension of the action space
dim = size(u_app,1);
% number of states
dim_x = size(x,1);
% set horizon size
horizon = 10;
% number of discretizations from 0 to t_final
N = length(t);
h = t(2) - t(1);
% Scaling matrix Sw
Sw = diag([1,1,1]);
STR.Sw = Sw;
% dimension of the context space : current state + desired state
dim_ctx = dim_x * 2; 
% number of samples in each example trajectory
num_samp = N-1;
% total number of samples in all trajectories
num_tot = num_samp * n;
% renaming t since it's used for iteration index in main loop
time = t;
% path is used only for IC
x_vec = x(:);
paths = zeros(dim_x * N, n);
% arrange u into proper GP input format
us = zeros(n,dim*N);
STR.h = h;
% adjust x by adding endpt to the horizon
x_hor = [x, repmat(x(:,end),1,horizon)];
% phi change limit
phi_lim = 0.2;

%% Get contexts and cost

if ~exist('ctx_TW.mat','file')
contexts = zeros(dim_ctx, num_samp, n);
cost = zeros(num_samp, n);
% to generate nominal prediction costs
handle = @robotTwoWheelsNominalKinematics;
% prepare us and costs
for i = 1:n
    % add ILC errors to get full path info
    paths(:,i) = x_vec + err(:,i);
    for j = 1:N-1
        % form the L vectors 
        u_lim = u_cnstr * ones(dim,1);
        x_cur = paths(dim_x*(j-1) + (1:dim_x),i);
        L = abs(h * robotTwoWheelsNominalKinematics(0,x_cur,PAR,u_lim,false));
        % get the optimized context from paths
        ctx_x = optim_trj(x_hor(1,j:j+horizon)', paths(dim_x*(j-1)+1,i), L(1));
        ctx_y = optim_trj(x_hor(2,j:j+horizon)', paths(dim_x*(j-1)+2,i), L(2));
        ctx_phi = optim_trj(x_hor(3,j:j+horizon)', paths(dim_x*(j-1)+3,i), phi_lim);
        ctx = [ctx_x'; ctx_y'; ctx_phi'];
        % convert into vector
        ctx_vec = ctx(:);
        % as context, we keep only present state + future state
        contexts(:,j,i) = ctx_vec(1:dim_ctx);
        % get the 2-norm difference
        dev = paths(dim_x*j+(1:dim_x),i) - ctx(:,2);
        cost(j,i) = dev'*(Sw'*Sw)*dev;
        % subtract nominal model prediction error
        next = step_RK4(h,paths(dim_x*(j-1)+(1:dim_x),i),u_app(:,j,i),CON,PAR,handle);
        % get nominal deviation
        pred_dev = next - ctx(:,2);
        pred_err = pred_dev'*(Sw'*Sw)*pred_dev;
        cost(j,i) = cost(j,i) - pred_err;
    end
end
save('ctx_TW.mat','contexts', 'cost');
else
load('ctx_TW.mat');
end

%% Fit hyperparameters on nominal model prediction error

%%%%%% TODO: DEBUG THIS! %%%%%%%%%%%
if ~exist('hyp_TW.mat','file')
% concat context space C with input space U
CxU = zeros(num_tot, dim_ctx + dim);
for i = 1:num_samp
    for j = 1:n
        idx = (i-1) *n + j;
        idx_u = dim * (i-1) + (1:dim);
        CxU(idx,:) = [contexts(:,i,j)', u_app(:,i,j)'];
    end
end
% format cost matrix into a vector
cost = cost(:,1:n);
cost = reshape(cost',num_tot,1);

%{ 
% set up a composite kernel (product of context, action kernels) 
maskCtx = [ones(1,dim_ctx), zeros(1,dim)]; % mask excluding action dimensions
maskAct = [zeros(1,dim_ctx), ones(1,dim)]; % mask excluding context dimensions
c1  = {@covLINard}; 
covCtx = {'covMask',{maskCtx,c1{:}}}; % context kernel
c2  = {@covSEiso};
covAct = {'covMask',{maskAct,c2{:}}}; % action kernel
covfunc = {@covProd,{covCtx,covAct}};  %covSum infers more noise
likfunc = @likGauss;
inf = @infExact;
meanfunc = []; %{@meanSum, {@meanLinear, @meanConst}}; 
%start with random hyperparameters (take exp for real values)
hyp.mean = []; %[zeros(dim+dim_ctx,1)];
hyp.cov = [rand(dim_ctx,1); 0; 0]; 
hyp.lik = log(0.5);
%fit hyperparameters
%function handle gp is used in training mode
%training: [nlZ dnlZ] = gp(hyp, inf, mean, cov, lik, x, y)
%nlZ: negative log probability of the training data
%dnlZ: derivative of above
Ncg = 100;
addpath('../gpml-matlab-v3.1-2010-09-27/util/');
hyp = minimize(hyp, @gp, -Ncg, inf, meanfunc, covfunc, likfunc, CxU, cost);
%inferred noise standard deviation
s = sprintf('Inferred noise standard deviation : %f', exp(hyp.lik));
disp(s);
%form them into struct so as to pass easily to function gpucb
GPSTR.hyp = hyp;
GPSTR.covfunc = covfunc;
GPSTR.likfunc = likfunc;
GPSTR.inf = inf;
GPSTR.meanfunc = meanfunc;
%}
% try different models
hyper; 
save('hyp_TW.mat', 'CxU', 'cost', 'GPSTR');
else
load('hyp_TW.mat');
end


%% Prepare variables for GP-UCB and plot no minal trajectory

% construct beta_t
% hoeffding prob. constant
delta = 0.1;
% these constants depend on the kernel bound (see Krause's GP-UCB paper)
a = 1; b = 1;
% bounds of the cube
r = u_cnstr - 1e-4;   
% dimension of the cube
d = dim_ctx + dim;
% beta multiplier - more aggression necessary for tackling real problems
beta_mult = 0.01;
beta_f = @(t,delta,d,a,b,r) beta_mult * (2*log((t^2)*2*(pi^2)/(3*delta)) +...
       2*d*log((t^2)*d*b*r*sqrt(log(4*d*a/delta))));
%iteration/run numbers
iter = N-1;
runs = 1; % each run learns from the previous runs
%regr = zeros(iter,1);
%cumr = zeros(iter,1);
% Establish bounds for variables
STR.bounds = r * [-ones(dim,1), ones(dim,1)];
% trajectory of the CGP-UCB
x_GP = zeros(dim_x, N, runs);
% real function 
fun_real = @robotTwoWheelsKinematics;
% nominal model
fun_nom = @robotTwoWheelsNominalKinematics;
STR.fun_nom = fun_nom;
% needed to calculate contexts (opt. trajectories)
u_lim = u_cnstr * ones(dim,1);

% plot desired trajectory and deviation if nominal input is applied
x_nom = sim_RK4(time,x,u_app(:,:,1),CON,PAR,fun_real);
for i = 1:size(x,1)
    subplot(2,3,i); 
    plot(time, x(i,:), '-r', time, x_nom(i,:), '-g');
    hold on;
end
subplot(2,3,[4 5 6]); 
plot(x(1,:), x(2,:), '-r', x_nom(1,:), x_nom(2,:), '-g');
hold on;

%% Bandit process
% maximize mu + sigma(x)

% if flag is set to 1, apply only nominal feedback
STR.FLAGS.nom_feedback = 1;
% for debugging
dbg_cost = zeros(1,iter);

for run = 1:runs
    
    % start from the initial point
    x_GP(:,1,run) = x(:,1);
    
    fprintf('Run number %d ...\n', run);
    %initialize cumulative regret vector
    %cumr1 = zeros(iter,1); 
    %sum1 = 0;
    STR.u_past = zeros(dim,1);
    STR.Kinv = 1;
            
    tic;
    for t = 1:iter
        % calculate the optimized context
        % form the L vectors 
        L = abs(h * robotTwoWheelsNominalKinematics(t,x_GP(:,t,run),PAR,u_lim,false));
        ctx_x = optim_trj(x_hor(1,t:t+horizon)', x_GP(1,t,run), L(1));
        ctx_y = optim_trj(x_hor(2,t:t+horizon)', x_GP(2,t,run), L(2));
        ctx_phi = optim_trj(x_hor(3,t:t+horizon)', x_GP(3,t,run), phi_lim);
        ctx = [ctx_x'; ctx_y'; ctx_phi'];
        ctx = ctx(:,[1,2]);
        STR.ctx = ctx(:);

        % beta is slightly increasing at each iteration
        STR.beta = beta_f((run-1) * iter + t,delta,d,a,b,r);
        
        % find new control signal to be tried
        % conditioning on the previous data is determined by past var.
        past = num_tot + (1: (run-1) * iter + t-1);
        STR.x_now = x_GP(:,t,run);
        u = gp_ucb_boost_R(GPSTR, CxU(past,:),cost(past),STR);
        %[u,STR] = gp_ucb_boost(t, CxU(past,:), cost(past), GPSTR.hyp, STR);
        % for fmincon to start from a nonzero value
        STR.u_past = u;

        % simulate one-step
        x_GP(:,t+1,run) = step_RK4(h,x_GP(:,t,run),u,CON,PAR,fun_real);
        dev = x_GP(:,t+1,run) - ctx(:,2);
        newCost = dev'*(Sw'*Sw)*dev;
        % subtract nominal model prediction error
        pred = step_RK4(h,x_GP(:,t,run),u,CON,PAR,fun_nom);
        pred_dev = pred - ctx(:,2);
        pred_err = pred_dev'*(Sw'*Sw)*pred_dev;
        newCost = newCost - pred_err;

        % add u and newCost to training data/observations pair
        if isempty(past), present = num_tot + 1;
        else present = past(end) + 1; %num_tot + t
        end
        CxU(present,:) = [ctx(:); u]'; 
        cost(present) = newCost; 
        
        % for debugging
        dF = x_GP(:,t+1,run) - pred;    
        dbg_cost(t) = dF'*(Sw'*Sw)*dF + 2*pred_dev'*(Sw'*Sw)*dF;
        
        % for regret calculation
        %sum1 = sum1 + cost(present);
        %cumr1(t) = sum1;
        
        % display every 10th 
        if(rem(t,10) == 0), fprintf('Time stage %d ...\n', t); end
    end
    
    STR.FLAGS.nom_feedback = 0;
    
    % display runtime
    fprintf('Trajectory simulation took %f seconds. \n', toc);
    
    % display total error in 2-norm
    err_total = norm(x_GP(1:2,:,run) - x(1:2,1:iter+1),2);
    fprintf('Total deviation is %f \n', err_total);
    
    % calculate average regret and avg. cum. regret
    %regr = ((run-1) * regr + cost(end-iter+1:end)) / run;
    %cumr = ((run-1) * cumr + cumr1) / run;
    
end

% plot followed trajectory and u
% plot final iteration result
for i = 1:size(x,1)
    subplot(2,3,i); 
    line = plot(time, x_GP(i,:,run));
end
subplot(2,3,[4 5 6]); plot(x_GP(1,:,run), x_GP(2,:,run), '-b');

% debug cost
figure;
plot(1:iter, cost(num_tot+1:end), '-r', 1:iter, dbg_cost);
legend('cost', 'debug');

% plot average regret
% figure(2)
% plot(1:iter, regr);
% title('(Instantaneous) Avg. Regret');
% 
% % plot average cumulative regret
% figure(3)
% plot(1:iter, cumr);
% title('Avg. Cumulative regret');