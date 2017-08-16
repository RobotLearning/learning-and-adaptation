%% diff. eq for reference tracking minimum principle with time varying A, B matrices

function xbar_dot = pmp_ltv_ref(t,xbar,params)

% xbar contains both x and lambda
n = length(xbar);
x = xbar(1:n/2);
lambda = xbar(n/2+1:end);

A = params.A;
B = params.B;
Q = params.Q;
R = params.R;
ts = params.ts;
ref = params.ref; % reference trajectory
r = interp1(ts,ref,t)';
x_dot = A*x - B*(R \ (B' * lambda));
lambda_dot = Q*(r-x) - A'*lambda;
xbar_dot = [x_dot(:); lambda_dot(:)];