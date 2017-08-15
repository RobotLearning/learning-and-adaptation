%% Testing batch inversion vs. recursive implementation

clc; clear; close all; rng(3);

% Testing the idea that minimum principle for trajectory tracking
% with quadratic tracking error costs and zero control input penalty
% gives the same results as batch lifted-matrix inversion

n = 2; % dim_x
m = 2; % dim_u
T = 1.0; % final time
N = 20; % num of traj. points
dt = T/N; % discretization
A = randn(n);
B = randn(n,m);
As = repmat(A,1,1,N);
Bs = repmat(B,1,1,N);
%A = randn(n,n,N);
%B = randn(n,m,N);

% Discretize and lift nominal  dynamics
[Ad,Bd] = discretizeDyn(As,Bs,dt);
F = liftDyn(Ad,Bd);

%% sample reference from a Gaussian Process

hp.type = 'squared exponential ard';
hp.l = 0.2;
hp.scale = 1;
hp.noise.var = 0.01;
t = dt * (1:N);
l = hp.l;
s = hp.scale;
mu = zeros(N,n);
kern = @(x1,x2) s * exp(-(norm(x1(:)-x2(:),2)^2)/(2*(l^2)));
Sigma = ker_matrix(t,kern);
[U,S] = eig(Sigma);
r = mu + U * sqrt(max(0,real(S))) * randn(N,n);

%plot(x,r);

%% Invert the reference to find best controls and evolve system from x0=0

rr = r';
r_lift = rr(:);
u_lift = pinv(F) * r_lift;
x_lift = F * u_lift;
x = reshape(x_lift,n,N);
plot(t,x,t,r);
disp('Max tracking error:');
max(max(x - r'))
disp('Due to inversion error:');
norm(eye(N*n) - F*pinv(F), 2)

%% Solve the BVP from MP

sol = bvp4c(odefun,bcfun,solinit)