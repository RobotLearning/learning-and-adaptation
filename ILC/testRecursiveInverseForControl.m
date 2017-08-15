%% Testing batch inversion vs. recursive implementation

clc; clear; close all; rng(3);

% Testing the idea that minimum principle for trajectory tracking
% with quadratic tracking error costs and small control input penalty
% gives the same results as batch lifted-matrix inversion as R \to 0

n = 2; % dim_x
m = 2; % dim_u
T = 1.0; % final time
N = 20; % num of traj. points
dt = T/N; % discretization
A = randn(n);
B = randn(n,m);
As = repmat(A,1,1,N+1);
Bs = repmat(B,1,1,N+1);
%A = randn(n,n,N);
%B = randn(n,m,N);

% Discretize and lift nominal dynamics
[Ad,Bd] = discretizeDyn(As,Bs,dt);

%% sample reference from a Gaussian Process

hp.type = 'squared exponential ard';
hp.l = 0.2;
hp.scale = 1;
hp.noise.var = 0.01;
t = dt * (1:N+1);
l = hp.l;
s = hp.scale;
mu = zeros(N+1,n);
kern = @(x1,x2) s * exp(-(norm(x1(:)-x2(:),2)^2)/(2*(l^2)));
Sigma = ker_matrix(t,kern);
[U,S] = eig(Sigma);
r = mu + U * sqrt(max(0,real(S))) * randn(N+1,n);

%plot(x,r);

%% Invert the reference to find best controls and evolve system from x0=r0
% TODO : check where traj starts!
rr = r';
r_lift = rr(:);
F = liftDyn(Ad,Bd);
u_lift = F \ r_lift;
x_lift = F * u_lift;
x = reshape(x_lift,n,N+1);
plot(t,x,t,r);
disp('Max tracking error:');
max(max(x - r'))
disp('Due to inversion error:');
norm(eye(N*n) - F*pinv(F), 2)

%% Solve the BVP from MP

sol = bvp4c(odefun,bcfun,solinit)
norm(eye((N+1)*n) - F*pinv(F), 2)

%% Solve u's instead with BVP (minimum principle)

solinit = bvpinit(linspace(0,T,N+1),zeros(1,2*n));

params.As = Ad;
params.Bs = Bd;
params.Q = eye(n);
params.R = 0.001*eye(m);
params.ref = r;
params.T = T;

odefun = @(t,x) pmp_ltv_ref(t,x,params);
bcfun = @(xa,xb) [xa(1:n)-rr(:,1); xb(n+1:end);];
sol = bvp4c(odefun,bcfun,solinit);
%u = -params.R \ (params.Bs' * ...
