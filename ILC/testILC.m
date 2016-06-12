%% Apply ILC

clc; clear; close all;
rng(3);

n = 4; % dim_x
m = 2; % dim_u
T = 1.0; % final time
N = 50;
dt = T/N; % discretization
A_nom = randn(n,n,N); % Nominal system
B_nom = randn(n,m,N);

% Discretize and lift nominal  dynamics
[Ad_nom,Bd_nom] = discretizeDyn(A_nom,B_nom,dt);
F_nom = liftDyn(Ad_nom,Bd_nom);

% Generate perturbation matrix and the actual dynamics
alpha = 10*min(svd(F_nom));
F_act = F_nom + genPerturbationMatrix(n,m,N,alpha);

checkAsymptoticStability(F_act,F_nom);

% observe errors
K = 100; % number of trials / experiments
E = zeros(N*n,K); % error along trajectory for each trial
U = zeros(N*m,K); % inputs along trajectory for each trial
err_norm = zeros(1,K);

% Reference trajectory is a draw from a Gaussian Process
traj = sampleTrajectory(n,N);
ref = traj';
ref = ref(:);

disp('Min error with actual model pinv');
err_min_norm = norm((eye(n*N)-F_act*pinv(F_act))*ref,2)
disp('Min error with nom model pinv');
err_min_nom_norm = norm((eye(n*N)-F_act*((pinv(F_nom)*F_act)\pinv(F_nom)))*ref,2)

% we assume no noise for now
for i = 1:K-1
    % get error
    E(:,i) = F_act * U(:,i) - ref;
    % ILC happening here
    U(:,i+1) = ilc(U(:,1:i),E(:,1:i),'simple',F_nom);
    % get norm of error
    err_norm(i) = norm(E(:,i),2);
end
% get last error
E(:,K) = F_act * U(:,K) - ref;
err_norm(K) = norm(E(:,K),2);

plot(err_norm);
ylabel('Error norm');
xlabel('Iterations');
title('Error norm vs trials');
