%% Maximum Likelihood Estimation of Hyperparameters
%
% Loads datasets generated using generation.m and then performs MLE on the
% data to estimate parameters (covariance with/without mean parameters). 
% Using Rasmussen's toolbox for gp and minimize functions.
%
% Author: Okan Koc, IDSC, ETH Zurich
%
% TODO: check REML!
%--------------------------------------------------------------------------


clc; clear; close all; diary('estimationResults.txt');

% define model structure
% mean not included
meanfunc = []; hyp.mean = [];
% covariance function hyperparameters: length ell and scale sf
covfunc = {@covSEiso};
% noise included
likfunc = @likGauss; 

% load datasets
load('estimation1D.txt');
load('hp.mat'); % exact parameter values
x = estimation1D(:,1);
% zero-mean values
y = estimation1D(:,2:6);
% mean values added
y2 = estimation1D(:,7:11);
% initial value for parameters
p0 = 1;
% initialize variables and structs
GPSTR.covfunc = covfunc;
GPSTR.likfunc = likfunc;
GPSTR.inf = @infExact;
GPSTR.meanfunc = meanfunc;
GPSTR.hyp = hyp;
ell1s = linspace(0.05,1,20);
ell2s = linspace(0.05,1,20);
sfs = linspace(0.4,1.2,9);
sns = linspace(0.1,0.5,5);
beta0s = linspace(-0.2,1,13);
beta1s = linspace(-0.1,1.5,17);
beta2s = linspace(-1.0,1.3,24);
nll3 = zeros(length(ell1s),length(sfs),length(sns));
% number of line searches for minimization
Ncg = 1000;


%% Sampling Based Optimization & Derivative Based Optimization

disp('Results for zero mean 1D');
disp('Covariance parameters maximizing the log likelihood:');

for k = 1:size(y,2)
    % To find optimal parameters
    % negative log likelihoods of different models
    fprintf('Dataset %d : \n', k);
    disp('Sampling based optimization:');
    
    % Using DIRECT 
    % Send options to Direct 
    options.showits   = 0;
    options.tol       = 0.05;
    options.maxevals  = 500;
    options.maxits    = 250;
    options.maxdeep   = 100;
    % Pass function as part of a Matlab Structure

    lb = [ell1s(1); sfs(1); sns(1)];
    ub = [ell1s(end); sfs(end); sns(end)];
    bounds = [lb, ub];
    problem.f = @(hyp) gps(GPSTR, x, y(:,k), hyp);
    [~,hpm,hist] = Direct(problem,bounds,options);
    fprintf('ell = %f, sf = %f, sn = %f \n', hpm(1), hpm(2), hpm(3));

    %{ 
    %plot the hyperparameter surface
    figure;
    surf(sfs,ells,-nlls);
    title('Log likelihood of the data under different parameters');
    xlabel('sf parameter values');
    ylabel('ell parameter values');
    %}     
    % Derivative Based Optimization
    % run conjugate gradient method to find best set of parameters

    %function handle gp is used in training mode
    %training: [nlZ dnlZ] = gp(hyp, inf, mean, cov, lik, x, y)
    %nlZ: negative log probability of the training data
    %dnlZ: derivative of above
    
    disp('MLE: Derivative based optimization...');
    fprintf('Starting from values ell1_0 = %f, sf0 = %f, sn0 = %f \n', ...
         p0, p0, p0);
    hyp.cov = log([p0; p0]);
    hyp.lik = log(p0);
    hyp = minimin(hyp, @gp, -Ncg, false, @infExact, meanfunc, covfunc, ...
                   likfunc, x, y(:,k));
    %inferred values
    fprintf('ell = %f, sf = %f, sn = %f \n', ...
             exp(hyp.cov(1)), exp(hyp.cov(2)), exp(hyp.lik));
    fprintf('True values: ell = %f, sf = %f, sn = %f \n', ...
             hp(k,1), hp(k,2), hp(k,3));
end


%% Adding Nonzero Mean
% complicates things!

ml = {@meanLinear};
mc = {@meanConst};
meanfunc2 = {'meanSum',{mc,ml}};
GPSTR.meanfunc = meanfunc2;

fprintf('\n');
disp('Results for 1D linear mean MLE:');
for k = 1:size(y2,2)
    % To find optimal parameters
    % negative log likelihoods of different models
    fprintf('Dataset %d : \n', k);
    disp('Sampling based optimization:');
    
    % Using DIRECT 
    % Send options to Direct 
    options.showits   = 0;
    options.tol       = 0.05;
    options.maxevals  = 2000;
    options.maxits    = 1000;
    options.maxdeep   = 100;
    % Pass function as part of a Matlab Structure

    lb = [ell1s(1); sfs(1); sns(1); beta0s(1); beta1s(1)];
    ub = [ell1s(end); sfs(end); sns(end); beta0s(end); beta1s(end)];
    bounds = [lb, ub];
    problem.f = @(hyp) gps(GPSTR, x, y2(:,k), hyp);
    [~,hpm,hist] = Direct(problem,bounds,options);
    fprintf('ell = %f, sf = %f, sn = %f \n', hpm(1), hpm(2), hpm(3));
    fprintf('mean slope = %f, intercept = %f \n', ...
             hpm(4), hpm(5));

    % Derivative Based Optimization
    % run conjugate gradient method to find best set of parameters
    disp('MLE: Derivative based optimization...');
    fprintf('Starting from values ell1_0 = %f, sf0 = %f, sn0 = %f \n', ...
         p0, p0, p0);
    fprintf('beta0 = %f, beta1 = %f \n', p0, p0);
    hyp.cov = log([p0; p0]);
    hyp.lik = log(p0);
    hyp.mean = [p0; p0];
    hyp = minimin(hyp, @gp, -Ncg, false, @infExact, meanfunc2, covfunc, ...
                   likfunc, x, y2(:,k));
    %inferred values
    fprintf('ell = %f, sf = %f, sn = %f \n', ...
             exp(hyp.cov(1)), exp(hyp.cov(2)), exp(hyp.lik));
    fprintf('mean slope = %f, intercept = %f \n', ...
             hyp.mean(2), hyp.mean(1));
         
    % REML !!!
    disp('REML...');
    fprintf('Starting from values ell1_0 = %f, sf0 = %f, sn0 = %f \n', ...
         p0, p0, p0);
    fprintf('beta0 = %f, beta1 = %f \n', p0, p0);
    hyp.cov = log([p0; p0]);
    hyp.lik = log(p0);
    hyp.mean = [p0; p0];
    [hyp,bhat] = reml(hyp, covfunc, likfunc, x, y2(:,k));
    %inferred values
    fprintf('ell = %f, sf = %f, sn = %f \n', ...
             exp(hyp.cov(1)), exp(hyp.cov(2)), exp(hyp.lik));
    fprintf('mean slope = %f, intercept = %f \n', ...
             bhat(2), bhat(1));
    fprintf('True values: ell = %f, sf = %f, sn = %f \n', ...
             hp(k,1), hp(k,2), hp(k,3));
    fprintf('mean slope = %f, intercept = %f \n', ...
             hp(k,5), hp(k,4));
    
end

%% Repeat in R2

% define model structure
% mean not included
meanfunc = [];
GPSTR.meanfunc = meanfunc; 
hyp.mean = [];
% covariance function hyperparameters: length ell and scale sf
covfunc = {@covSEard};
GPSTR.covfunc = covfunc;

% load datasets
load('estimation2D.txt');
load('hp2.mat'); % exact parameter values
X = estimation2D(:,1:2);
% zero-mean values0
y = estimation2D(:,3:7);
% mean values added
y2 = estimation2D(:,8:12);

fprintf('\n');
disp('Results for zero mean 2D');
disp('Covariance parameters maximizing the log likelihood:');

for k = 1:size(y,2)
    % To find optimal parameters
    % negative log likelihoods of different models
    fprintf('Dataset %d : \n', k);
    disp('Sampling based optimization:');
    
    % Using DIRECT 
    % Send options to Direct 
    options.showits   = 0;
    options.tol       = 0.05;
    options.maxevals  = 1000;
    options.maxits    = 500;
    options.maxdeep   = 100;
    % Pass function as part of a Matlab Structure

    lb = [ell1s(1); ell2s(1); sfs(1); sns(1)];
    ub = [ell1s(end); ell2s(end); sfs(end); sns(end)];
    bounds = [lb, ub];
    problem.f = @(hyp) gps(GPSTR, X, y(:,k), hyp);
    [~,hpm,hist] = Direct(problem,bounds,options);
    fprintf('ell1 = %f, ell2 = %f, sf = %f, sn = %f \n', ...
             hpm(1), hpm(2), hpm(3), hpm(4));

    % Derivative Based Optimization
    % run conjugate gradient method to find best set of parameters    
    disp('MLE: Derivative based optimization...');
    fprintf('Starting from values ell1_0 = %f, ell2_0 = %f, sf0 = %f, sn0 = %f \n', ...
         p0, p0, p0, p0);
    hyp.cov = log([p0; p0; p0]);
    hyp.lik = log(p0);
    hyp = minimin(hyp, @gp, -Ncg, false, @infExact, meanfunc, covfunc, ...
                   likfunc, X, y(:,k));
    %inferred values
    fprintf('ell1 = %f, ell2 = %f, sf = %f, sn = %f \n', ...
             exp(hyp.cov(1)), exp(hyp.cov(2)), exp(hyp.cov(3)), ...
             exp(hyp.lik));
    fprintf('True values: ell1 = %f, ell2 = %f, sf = %f, sn = %f \n', ...
             hp2(k,1), hp2(k,2), hp2(k,3), hp2(k,4));
end

%% Including Mean in R2 case
% linear mean with intercept included
meanfunc2 = {@meanSum, {@meanConst, @meanLinear}}; 
GPSTR.meanfunc = meanfunc2;

fprintf('\n');
disp('Results for linear mean 2D MLE');
for k = 1:size(y2,2)
    % To find optimal parameters
    % negative log likelihoods of different models
    fprintf('Dataset %d : \n', k);
    disp('Sampling based optimization:');
    
    % Using DIRECT 
    % Send options to Direct 
    options.showits   = 0;
    options.tol       = 0.05;
    options.maxevals  = 8000;
    options.maxits    = 4000;
    options.maxdeep   = 100;
    % Pass function as part of a Matlab Structure

    lb = [ell1s(1); ell2s(1); sfs(1); sns(1); ...
          beta0s(1); beta1s(1); beta2s(1)];
    ub = [ell1s(end); ell2s(1); sfs(end); sns(end); 
          beta0s(end); beta1s(end); beta2s(1)];
    bounds = [lb, ub];
    problem.f = @(hyp) gps(GPSTR, X, y2(:,k), hyp);
    [~,hpm,hist] = Direct(problem,bounds,options);
    fprintf('ell1 = %f, ell2 = %f, sf = %f, sn = %f \n', ...
             hpm(1), hpm(2), hpm(3), hpm(4));
    fprintf('mean slope1 = %f, slope2 = %f, intercept = %f \n', ...
             hpm(6), hpm(7), hpm(5));

    % MLE: Derivative Based Optimization
    % run conjugate gradient method to find best set of parameters    
    disp('MLE: Derivative based optimization...');
    fprintf('Starting from values ell1_0 = %f, ell2_0 = %f, sf0 = %f, sn0 = %f \n', ...
         p0, p0, p0, p0);
    fprintf('beta0 = %f, beta1 = %f, beta2 = %f \n', ...
         p0, p0, p0);
    hyp.cov = log([p0; p0; p0]);
    hyp.lik = log(p0);
    hyp.mean = [p0; p0; p0];
    hyp = minimin(hyp, @gp, -Ncg, false, @infExact, meanfunc2, covfunc, ...
                   likfunc, X, y2(:,k));
    %inferred values
    fprintf('ell1 = %f, ell2 = %f, sf = %f, sn = %f \n', ...
             exp(hyp.cov(1)), exp(hyp.cov(2)), exp(hyp.cov(3)), ...
             exp(hyp.lik));
    fprintf('mean slope1 = %f, slope2 = %f, intercept = %f \n', ...
             hyp.mean(2), hyp.mean(3), hyp.mean(1));
         
    % REML !!!
    disp('REML...');
    fprintf('Starting from values ell1_0 = %f, ell2_0 = %f, sf0 = %f, sn0 = %f \n', ...
         p0, p0, p0, p0);
    fprintf('beta0 = %f, beta1 = %f, beta2 = %f \n', p0, p0, p0);
    hyp.cov = log([p0; p0; p0]);
    hyp.lik = log(p0);
    hyp.mean = [p0; p0; p0];
    [hyp,bhat] = reml(hyp, covfunc, likfunc, X, y2(:,k));
    %inferred values
    fprintf('ell1 = %f, ell2 = %f, sf = %f, sn = %f \n', ...
             exp(hyp.cov(1)), exp(hyp.cov(2)), exp(hyp.cov(3)), exp(hyp.lik));
    fprintf('mean slope1 = %f, slope2 = %f, intercept = %f \n', ...
             bhat(2), bhat(3), bhat(1));
    fprintf('True values: ell1 = %f, ell2 = %f, sf = %f, sn = %f \n', ...
             hp2(k,1), hp2(k,2), hp2(k,3), hp2(k,4));
    fprintf('mean slope1 = %f, slope2 = %f, intercept = %f \n', ...
             hp2(k,6), hp2(k,7), hp2(k,5));
end

diary off;