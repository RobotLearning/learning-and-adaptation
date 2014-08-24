%% HYPERPARAMETER ESTIMATION
%--------------------------------------------------------------------------
% Called once by the GP main script (GP_Main)
%
% Author: Okan Koc, IDSC Lab, ETH Zurich
%
% -------------------------------------------------------------------------

matfile = 'hyp_aircraft_deltaG.mat';
%matfile = 'hyp_aircraft_wind.mat';
savefile = [aircraftFolder, '/', matfile];
if ~exist(savefile,'file')

%%%%%%%%%%%%%%%%%%%%%%%%%% CONCAT CONTEXT WITH INPUT %%%%%%%%%%%%%%%%%%%%%%
% concat context space C with input space U

CxU = zeros(num_tot, dim_ctx + dim);
cost = zeros(num_tot, 1);
idx = 1;
for i = 1:n
    for j = 1:N-1
        CxU(idx,:) = [contexts{i}(:,j); us{i}(:,j)]';
        cost(idx) = costs{i}(j);
        idx = idx + 1;
    end
end

%%%%%%%%%%%%%%% RASMUSSENS TOOLBOX HYPERPARAMETER ESTIMATION %%%%%%%%%%%%%%
%%{
% set up a composite kernel (product of context, action kernels) 
%maskCtx = [0,0,0,1,1,0,0,0,0,1,0,0,zeros(1,dim)];
maskCtx = [ones(1,dim_ctx), zeros(1,dim)]; 
maskAct = [zeros(1,dim_ctx), ones(1,dim)];
c1 = {@covSEard};
c2  = {@covLINard};
covCtx = {'covMask',{maskCtx,c1{:}}}; % context kernel
covAct  = {'covMask',{maskAct,c2{:}}}; % action kernel
covfunc = {@covProd,{covCtx,covAct}};
%covfunc = {@covSum, {covCtx,covAct}};
likfunc = @likGauss;
inf = @infExact;
meanfunc = [];
hyp.mean = [];
hyp.cov = zeros(dim_ctx+dim+1,1);
%hyp.cov = zeros(3+dim+1,1);
hyp.lik = log(0.1);
%fit hyperparameters
%function handle gp is used in training mode
%training: [nlZ dnlZ] = gp(hyp, inf, mean, cov, lik, x, y)
%nlZ: negative log probability of the training data
%dnlZ: derivative of above
Ncg = 100;
addpath('../gpml-matlab-v3.1-2010-09-27/util/');
hyp = minimize(hyp, @gp, -Ncg, inf, meanfunc, covfunc, likfunc, CxU, cost);
%inferred noise standard deviation
fprintf('Inferred noise standard deviation : %f \n', exp(hyp.lik));
%form them into struct so as to pass easily to function gpucb
GPSTR.hyp = hyp;
GPSTR.covfunc = covfunc;
GPSTR.likfunc = likfunc;
GPSTR.inf = inf;
GPSTR.meanfunc = meanfunc;
%}
% try different models
%hyper; 
save(savefile, 'CxU', 'cost', 'GPSTR');

else
load(savefile);
end