%% Recursively estimate LTV system matrices with 
% Linear Bayesian Regression (LBR)

% assumption is that the trials are exactly the same
% for linear time varying system F
% and random inputs

clc; clear; close all;
m = 2;
n = 4;
N = 5;
K = N*m; % effective number of trials to estimate well = N*m
B = rand(n,m,N);
A = rand(n,n,N);
F = liftDyn(A,B);
u = rand(N*m,K);
e = F*u;

errNormRecursive = zeros(1,K);
errNormBatch = zeros(1,K);
errInvNorm = zeros(1,K);

% Generate perturbation matrix and the actual dynamics
alpha = 10*min(svd(F));
Fest = F + genPerturbationMatrix(n,m,N,alpha);
numNonzeroElem = sum(sum(F > 0));

% proportional to number of samples that was invested in 
% acquiring Fest (e.g. during system identification)
precision_val = 0.0001;
% precision matrix
Gamma = precision_val * eye(numNonzeroElem);
% prior inv gamma distribution values
hp.a = 2; 
hp.b = 0;

% for each trial
for i = 1:K
    [Fest,Gamma,hp] = recursiveEstSysMat(Fest,Gamma,hp,u(:,i),e(:,i),N);
    errNormRecursive(i) = norm(F-Fest);
    %errInvNorm(i) = norm(pinv(F)-pinv(Fest));
end

% check with batch estimation
for i = 1:K
    Fest2 = batchEstSysMat(u(:,1:i),e(:,1:i),[n,m]);
    errNormBatch(i) = norm(F-Fest2);
end

plot(1:K,errNormRecursive,1:K,errNormBatch);
legend('recursive','batch');
