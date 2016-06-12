%% Estimate LTV system matrices with least squares

% assumption is that the trials are exactly the same
% for linear time varying system F
% and random inputs

clc; clear; close all;
m = 2;
n = 2;
N = 2;
K = N*m; % effective number of trials to estimate well = N*m
B = rand(n,m,N);
A = rand(n,n,N);
F = liftDyn(A,B);
u = rand(N*m,K);
e = F*u;

tic
systemSize = [n,m];
Fest = estimateSystemMatrices(u,e,systemSize);
toc

errNorm = norm(F-Fest)
errInvNorm = norm(pinv(F)-pinv(Fest))

% try something else
% create the duplication matrix
% tic
% T = tril(ones(N));
% T = kron(T(:)',ones(N));
% T = [T(1,1:4),T(2,1:4), T(1,5:end), T(2,5:end)];
% TT = diag(T);
% D2 = TT(logical(T),:)';
% M2 = kron(u',eye(n*N))*D2;
% E2 = e(:);
% Est2 = pinv(M2)*E2;
% Est2AddZero = D2*Est2;
% F_est2 = reshape(Est2AddZero',n*N,m*N);
% toc