%% HYPERPARAMETER ESTIMATION
%--------------------------------------------------------------------------
% Called once by the GP main script (GP_Main)
%
% Author: Okan Koc, IDSC Lab, ETH Zurich
%
% -------------------------------------------------------------------------

%matfile = 'hyp_quad_actuator.mat';
matfile = 'hyp_quad_wind.mat';
savefile = strcat(qcopterFolder, '/', matfile);
if ~exist(savefile,'file')

%%%%%%%%%%%%%%%%%%%%%%%%%% CONCAT CONTEXT WITH INPUT %%%%%%%%%%%%%%%%%%%%%%
% concat context space C with input space U
CxU = zeros(num_tot, dim_ctx + dim);
cost = [];
idx = 1;
for i = 1:n
    % format cost cell into a vector
    cost = [cost; costs{i}]; %#ok
    for j = 1:N(i)-1
        CxU(idx,:) = [contexts{i}(:,j); us{i}(:,j)]';
        idx = idx + 1;
    end
end

%%%%%%%%%%%%%%% RASMUSSENS TOOLBOX HYPERPARAMETER ESTIMATION %%%%%%%%%%%%%%
%%{
% set up a composite kernel (product of context, action kernels) 
maskCtx = [1,1,0,0,0,1,1,1,1,0, zeros(1,dim)]; % wind case
%maskCtx = [0,1,0,0,1,0,1,0,0,0,zeros(1,dim)]; %actuator case
maskAct = [zeros(1,dim_ctx), 1, 0];
c1 = {@covSEard};
covCtx = {'covMask',{maskCtx,c1{:}}}; % context kernel
c2  = {@covLIN};
%c2 = {'covPoly', 2};
covAct = {'covMask',{maskAct,c2{:}}}; % action kernel
covfunc = {@covProd,{covCtx,covAct}};  %covSum infers more noise

likfunc = @likGauss;
inf = @infExact;
meanfunc = []; 
%start with random hyperparameters (take exp for real values)
hyp.mean = [];
hyp.cov = [zeros(sum(maskCtx),1);log(0.0064)];
%hyp.cov = [zeros(sum(maskCtx),1);log(0.0064);0;0];
hyp.lik = log(0.001);
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
save(savefile, 'CxU', 'cost', 'GPSTR');
%}

else
load(savefile);
end