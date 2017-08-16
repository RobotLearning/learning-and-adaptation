%% Testing batch inversion vs. recursive implementation

clc; clear; close all; rng(3);

% Testing the idea that minimum principle for trajectory tracking
% with quadratic tracking error costs and small control input penalty
% gives the same results as batch lifted-matrix inversion as R \to 0

n = 3; % dim_x
m = 1; % dim_u
T = 1.0; % final time
N = 50; % num of traj. points
dt = T/N; % discretization
A = randn(n); %[0 1 0; 0 0 1; -4 -6 -4];
B = randn(n,m); %[0; 0; 1];
As = repmat(A,1,1,N);
Bs = repmat(B,1,1,N);
%A = randn(n,n,N);
%B = randn(n,m,N);

%% sample reference from a Gaussian Process

hp.type = 'squared exponential ard';
hp.l = 0.2;
hp.scale = 1;
hp.noise.var = 0.01;
t = dt * (0:1:N);
l = hp.l;
s = hp.scale;
mu = zeros(N+1,n);
kern = @(x1,x2) s * exp(-(norm(x1(:)-x2(:),2)^2)/(2*(l^2)));
Sigma = ker_matrix(t,kern);
[U,S] = eig(Sigma);
r = mu + U * sqrt(max(0,real(S))) * randn(N+1,n);
r = r - repmat(r(1,:),N+1,1); % make sure traj starts from 0
t = t(2:end);
r = r(2:end,:);

%plot(x,r);

%% Invert the reference to find best controls and evolve system from x0=0

rr = r';
r_lift = rr(:);
% Discretize and lift nominal dynamics
tic; 
[Ad,Bd] = discretizeDyn(As,Bs,dt); 
toc
F = liftDyn(Ad,Bd); % seems this operation is costly
toc
tic
u_lift = F \ r_lift;
fprintf('Plant inversion took %f sec.\n', toc);
x_lift = F * u_lift;
u_invert = reshape(u_lift,m,N);
x_invert = reshape(x_lift,n,N);
figure(1);
plot(t,x_invert,t,r);
disp('Max tracking error:');
max(max(x_invert - r'))
disp('Due to inversion error:');
norm(eye(N*n) - F*pinv(F), 2)
disp('Max control');
max(max(abs(u_invert)))

%% Solve u's instead with BVP (minimum principle)

tic;
solinit = bvpinit(linspace(dt,T,N),zeros(1,2*n));
params.A = A;
params.B = B;
params.Q = eye(n);
params.R = 0.0001*eye(m);
params.ref = r;
params.ts = t;
odefun = @(t,x) pmp_ltv_ref(t,x,params);
bcfun = @(xa,xb) [xa(1:n)-rr(:,1); xb(n+1:end);];
sol = bvp4c(odefun,bcfun,solinit);
sol_at_t = deval(sol,t);
x_incr = sol_at_t(1:n,:);
lambda = sol_at_t(n+1:end,:); % TODO: should be evaluated from 0 to T-dt
u_incr = -params.R\(B'*lambda);
fprintf('Stable (incr.) inversion took %f sec.\n', toc);
figure(2);
plot(t,x_incr,t,r);
disp('Max tracking error:');
max(max(x_incr - r'))
disp('Max control');
max(max(abs(u_incr)))