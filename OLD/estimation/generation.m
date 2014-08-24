%% Generate zero-mean functions drawn from GPs 
% with different hyperparameters
% generate a syntetic data set drawn from a known exponential kernel

clear all, close all; clc;

% define model structure
% hyperparameter sets are defined here
% ell, sf and sn, intercept and slope respectively
hp = [0.150, 0.5, 0.1, 1, 0.5;
      0.250, 0.8, 0.2, 0.4, 0.2;
      0.100, 1.0, 0.3, 0.6, 0.1;
      0.300, 0.5, 0.1, 0.8, 1.0;
      0.150, 1.2, 0.2, 0.2, 1.5];
numdataset = size(hp,1);
% mean not included
meanfunc = []; hyp.mean = [];
% covariance function hyperparameters: length ell and scale sf
covfunc = {@covSEiso};
% noise included
likfunc = @likGauss; 

% number of input points
n = 50;
%rng('default');
x = linspace(0,1,n); x = x(:);
%x = sort(x(:)); % useful if you randomize input points

%% Generate observations for each set

y = zeros(n, numdataset);
% add means to y2
y2 = zeros(n, numdataset);

for i = 1:numdataset
    ell = hp(i,1); sf = hp(i,2);
    hyp.cov = log([ell; sf]);
    sn = hp(i,3);
    hyp.lik = log(sn);
    beta1 = hp(i,5);
    beta0 = hp(i,4);

    % covariance matrix 
    K = feval(covfunc{:}, hyp.cov, x);
    % mean values 
    mu = beta1 * x + beta0;
    % statistical toolbox function mvnrnd uses eigendecomposition - more stable
    %y = mvnrnd(mu, K); y = y(:);
    %y = y + exp(hyp.lik)*randn(n,1); % add noise
    % usual method for generating y - 
    % for f it does not work for large n unless noise is added
    Kn = K + (exp(hyp.lik)^2)*eye(length(K));
    y(:,i) = chol(Kn)'*randn(n, 1);
    y2(:,i) = y(:,i) + mu;
    
    % plot the function
    figure(i);
    %s2 = ['Zero mean fnc values drawn from a sq-exp kernel with par: ', s];
    plot(x, y, '+');
    %xlabel('x');
    %ylabel('y values');
    %title(s2);

    % draw the mean and variances of test points
    z = linspace(min(x), max(x), 100)';
    [m, s2] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, y(:,i), z);
    f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)]; 
    fill([z; flipdim(z,1)], f, [7 7 7]/8)
    hold on; 
    plot(z, m, '-r'); 
    plot(x, y(:,i), '+');
end

%% save functions in a R-compatible ASCII data file

xy = [x, y, y2];
save('estimation1D.txt', 'xy', '-ASCII', '-tabs');
save('hp.mat', 'hp');


%% Increasing Dimensionality of the Function
% Automatic Relevance Determination (ARD) in R2

clc; clear all; close all;
% define model structure
% hyperparameter sets are defined here
% ell1, ell2, sf, sn, intercept, slope1 and slope2 respectively
hp2 = [0.150, 0.5, 0.5, 0.1, 1.0, 0.5, -1.0;
       0.250, 0.35, 0.4, 0.2, 0.4, 0.8, 0.2;
       0.100, 0.8, 1.0, 0.3, 1.0, -0.1, 0.6;
       0.300, 0.15, 0.5, 0.1, -0.5, 1.0, 0.8;
       0.150, 0.2, 1.2, 0.1, -0.2, 1.5, 1.3];
numdataset2 = size(hp2,1);
% mean not included
meanfunc = []; hyp.mean = [];
% covariance function
covfunc = {@covSEard}; 
% noise included
likfunc = @likGauss; 

% generate new function values corrupted with noise
n = 10;
%rng('default');
x = linspace(0,1,n);
% generate n x n training points
xx = repmat(x,n,1);
xxx = xx(:);
yy = xx';
yyy = yy(:);
X = [xxx, yyy]; % design matrix / training points

%% Generate observations for each set

y = zeros(n^2, numdataset2);
% add means to y2
y2 = zeros(n^2, numdataset2);

for i = 1:numdataset2
    ell = hp2(i,1:2); sf = hp2(i,3);
    hyp.cov = log([ell'; sf]);
    sn = hp2(i,4);
    hyp.lik = log(sn);
    beta = hp2(i,6:7);
    beta0 = hp2(i,5);

    % covariance matrix 
    K = feval(covfunc{:}, hyp.cov, X);
    % mean values 
    mu = X * beta(:) + beta0;
    % statistical toolbox function mvnrnd uses eigendecomposition - more stable
    %y = mvnrnd(mu, K); y = y(:);
    %y = y + exp(hyp.lik)*randn(n,1); % add noise
    % usual method for generating f - numerically problematic for large n
    Kn = K + (exp(hyp.lik)^2)*eye(length(K));
    y(:,i) = chol(Kn)'*randn(n*n, 1);
    y2(:,i) = y(:,i) + mu;
    
    % plot the function
    figure(i);
    [Xs,Ys] = meshgrid(x,x);
    Z = reshape(y(:,i), n, n);
    surf(Xs,Ys,Z);
    title('Function values drawn from a kernel with known parameters');
end

% save functions in a R-compatible ASCII data file

Xy = [X, y, y2];
save('estimation2D.txt', 'Xy', '-ASCII', '-tabs');
save('hp2.mat', 'hp2');