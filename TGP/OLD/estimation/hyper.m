%% Hyperparameter estimation with Rasmussen's toolbox (GP)

% This script picks among different models one with the best likelihood
Ncg = 100;
addpath('../gpml-matlab-v3.1-2010-09-27/util/');

% ways of forming kernel covariances
covfuncs = {@covSEiso, @covLIN, @covSEard, @covLINard};
param_covfuncs = {rand(2,1), [], rand(dim_ctx+1,1), rand(dim_ctx,1)};
% ways of forming mean functions
meanfuncs = {@meanZero, @meanConst, @meanLinear};
param_meanfuncs = {[], rand, rand(dim+dim_ctx,1)};
% ways of forming composite covariances
covcomps = {@covSum, @covProd};

% set up a composite kernel (product of context, action kernels) 
maskCtx = [ones(1,dim_ctx), zeros(1,dim)]; % mask excluding action dimensions
maskAct = [zeros(1,dim_ctx), ones(1,dim)]; % mask excluding context dimensions

% we want to minimize the negative log likelihood
nlZ = zeros(1,16);
% index all the methods to be tried
i = 1;
for covsctx = 3:4 % covariance candidates for context space
    for covsact = 1:2 % action space covariances
        for means = 1:2 % mean functions to be tried
            for comp = 1:2 % ways of forming composite covariances
                c1 = covfuncs{covsctx};
                covCtx = {'covMask',{maskCtx,c1}}; % context kernel
                c2 = covfuncs{covsact};
                covAct = {'covMask',{maskAct,c2}}; % action kernel
                covfunc = {covcomps{comp},{covCtx,covAct}};  %covSum infers more noise
                likfunc = @likGauss;
                inf = @infExact;
                meanfunc = meanfuncs{means};
                %start with random hyperparameters (take exp for real values)
                hyp.mean = param_meanfuncs{means};
                hyp.cov = [param_covfuncs{covsctx}; param_covfuncs{covsact}];
                hyp.lik = log(0.1);
                hyp = minimize(hyp, @gp, -Ncg, ...
                               inf, meanfunc, covfunc, likfunc, CxU, cost);
                %form them into struct so as to pass easily to function gpucb
                GPSTR(i).hyp = hyp;
                GPSTR(i).covfunc = covfunc;
                GPSTR(i).likfunc = likfunc;
                GPSTR(i).inf = inf;
                GPSTR(i).meanfunc = meanfunc;
                s = sprintf('Inferred noise standard deviation : %f', exp(hyp.lik));
                disp(s);

                [nlZ(i) dnlZ] = gp(hyp, inf, meanfunc, covfunc, likfunc, CxU, cost);
                i = i+1;
            end
        end
    end
end

[~,best] = min(nlZ);
GPSTR = GPSTR(best);