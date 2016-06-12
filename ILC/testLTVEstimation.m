%% Estimate LTV system matrices with least squares

% assumption is that the trials are exactly the same
% for linear time varying system F
% and random inputs

clc; clear; close all;
m = 2;
n = 3;
N = 10;
rng(2)
K = N*m; % effective number of trials to estimate well = N*m
B = rand(n,m,N);
A = rand(n,n,N);
F = liftDyn(A,B);
u = rand(N*m,K);
e = F*u;

tic
systemSize = [n,m];
Fest = batchEstSysMat(u,e,systemSize);
toc

errNorm = norm(F-Fest)
% errInvNorm = norm(pinv(F)-pinv(Fest))

tic 
% another method that does not have any for loops 
%but doesnt consider block structure
[D,E] = genDuplicationMatrix(N,n,m);
M = kron(u',eye(n*N))*D;
Est = pinv(M)*e(:);
EstAddZero = D*Est;
Fest2 = reshape(EstAddZero',n*N,m*N);
toc

errNorm = norm(F-Fest2)