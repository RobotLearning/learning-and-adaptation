%% diff. eq for reference tracking minimum principle with time varying A, B matrices

function xbar_dot = pmp_ltv_ref(t,xbar,params)

% xbar contains both x and lambda
n = length(xbar);
x = xbar(1:n/2);
lambda = xbar(n/2+1:end);

As = params.As;
Bs = params.Bs;
Q = params.Q;
R = params.R;
ref = params.ref'; % reference trajectory
T = params.T; % final time tf
N = size(As,3)-1;
dt = T/N;
idx_floor = min(N,floor(t/dt)+1);
t_rem = t - floor(t/dt)*dt;


A = As(:,:,idx_floor);
B = Bs(:,:,idx_floor);
r = ref(:,idx_floor);
if t_rem > 0.0
    A = A + (t_rem/dt) * As(:,:,idx_floor+1); 
    B = B + (t_rem/dt) * Bs(:,:,idx_floor+1); 
    r = r + (t_rem/dt) * ref(:,idx_floor+1);
end

x_dot = A*x - R \ (B' * lambda);
lambda_dot = Q*(r-x) - A'*lambda;
xbar_dot = [x_dot(:); lambda_dot(:)];