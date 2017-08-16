%% Solving Optimal Control with Gradient descent 

% Randomize controls
% Evolve x with those controls 
% Evolve lambda back with those x
% Now use derivative of Hamiltonian to descend a step

n = 1; % dim_x
m = 1; % dim_u
T = 1.0; % final time
N = 50; % num of traj. points
dt = T/N; % discretization
A = randn(n); %[0 1 0; 0 0 1; -4 -6 -4];
B = randn(n,m); %[0; 0; 1];
Q = eye(n);
R = 0.001*eye(m);
QQ = []; RR = [];
for l = 1:N
    QQ = blkdiag(QQ,Q);
    RR = blkdiag(RR,R);
end

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
r = r';

%% Solve with Gradient Descent
u = zeros(m,N);
x0 = zeros(n,1);
x = zeros(n,N+1);
lambda = zeros(n,N+1);
x(:,1) = x0;
lambda(:,end) = 0;
alpha = zeros(1,N);
num_iter = 20;
for i = 1:num_iter
    
    % evolve x
    for j = 1:N
        x(:,j+1) = A*x(:,j) + B*u(:,j);
    end
    % evolve back lambda
    for j = N:-1:1
        lambda(:,j) = Q*(x(:,j+1) - r(:,j+1)) + A'*lambda(:,j+1);
    end    
    dHdu = R*u + B'*lambda(:,1:end-1);
    % take a step
    for k = 2:N+1
        f = @(u) deal(0.5*(x(:,k)-r(:,k))'*Q*(x(:,k)-r(:,k)) + ...
                 0.5*u'*R*u + lambda(:,k)'*(A*x(:,k)+B*u),...
                      R*u + B'*lambda(:,k)); % grad
        alpha(k-1) = linesearch(f,-dHdu(:,k-1),u(:,k-1),0.5,10^(-4));
        u(:,k-1) = u(:,k-1) - alpha(k-1)*dHdu(:,k-1);
    end
    % u = u - 0.001*dHdu(:)'; % BLOWS UP!
    % give info on total cost
    J(i) = 0.5*(x(2:end)-r(2:end))*QQ*(x(2:end)-r(2:end))' + 0.5*u*RR*u';
end
J