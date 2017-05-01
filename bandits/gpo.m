%
% Gaussian Process Optimization
%

function [idx,u,Kinv,s2max] = gpo(rewards,xs,Kinv,kern,mesh,acq)

t = size(rewards,2);
% if at start point return one of the endpoints
 if t == 0, 
     u = 0;
     idx = 1;
     s2max = 1;
     return; 
 end

% Update Inverse Matrix
% perform update for t
% calculate the matrix inverse incrementally
if t > 1 % if t is 1, then Kinv is a scalar
    % calculate A and store, calculate s
    Q = ker_vector(t-1,xs(t),xs(1:t-1),kern);
    A = Kinv * Q;
    S = kern(xs(t), xs(t)); % usually 1
    % calculate M
    m = 1/(S - Q'*A);
    % calculate Phat
    P_hat = Kinv + m*(A*A');
    % calculate Qhat
    Q_hat = -m*A;
    % form the full inverse matrix
    Kinv = [P_hat, Q_hat; Q_hat', m];
end

% Optimize the acquisition function
%{
options = optimset('Display','off','Algorithm','interior-point');
options = optimset(options,'GradObj','off');
% fmincon has the advantage of starting from a nonzero value
u0 = 1;
f = @(x) -acquisition_function(t,x,xs,Kinv,kern,rewards,acq);
u = fmincon(f,u0,[],[],[],[],lb,ub,[],options);
%}

% Taking maximum of the mesh
%%{

vals = acquisition_function_vec(t,mesh,xs,Kinv,kern,rewards,acq);
[~,idx] = max(vals);
u = mesh(idx);
s2max = calc_sigma_at_arm(t,u,xs,Kinv,kern);
%}

end

% Construct acquisition function
function val = acquisition_function(t,x,xs,invmat,kern,rewards,acq)

% calculate kernel vector k and store
ker_vec = ker_vector(t,x,xs,kern);
% create mu function 
mu = @(x) ker_vec' * (invmat * rewards(:));
%mu = @(x) ker_mat' * invmat * cost(:);

% create sigma function
% this makes sure we don't have negative values due to small num. errors
sigma2 = @(x) max(kern(x,x) - ker_vec' * (invmat * ker_vec),0);

% pass mu and sigma functions as function handles
val = acq(mu(x),sigma2(x));

end

function s2max = calc_sigma_at_arm(t,x,xs,invmat,kern)

% calculate kernel vector k and store
ker_vec = ker_vector(t,x,xs,kern);

% create sigma function
% this makes sure we don't have negative values due to small num. errors
sigma2 = @(x) max(kern(x,x) - ker_vec' * (invmat * ker_vec),0);

s2max = sigma2(x);

end

% Construct acquisition function over the whole mesh
function vals = acquisition_function_vec(t,x,xs,invmat,kern,rewards,acq)

% calculate kernel vector k and store
ker_mat = zeros(size(xs,2), size(x,2));
for i = 1:size(x,2)
    ker_mat(:,i) = ker_vector(t,x(:,i),xs,kern)';
end    
% calculate prior covariance
ker_prior_mat = ker_matrix(x,kern);
% create mu function 
mu = @(x) ker_mat' * (invmat * rewards(:));

% create sigma function
% this makes sure we don't have negative values due to small num. errors
cov = @(x) max(ker_prior_mat - ker_mat' * (invmat * ker_mat),0);

% pass mu and sigma functions as function handles
m = mu(x);
s2 = diag(cov(x));
vals = acq(m,s2);

end

%
% construct the kernel vector - covariances of the test pt and previous pts
% INPUTS: 
% t - last time stage
% x - test point
% us - points to be tried (up to horizon)
% fun - function handle for the kernel
% OUTPUT
% out - kernel column vector at the previous points and current test pt.
%
function out = ker_vector(t,x,us,fun)

Kvect = zeros(t,1);
for i = 1:t
    Kvect(i) = fun(us(:,i),x(:));
end
out = Kvect;
end

function out = ker_matrix(us,kernel)

%
% Constructs the kernel matrix all at once.
% Only used to generate test functions.
% INPUTS: 
% us - sampling points 
% kernel - handle for the kernel generating function
% OUTPUTS
% out - matrix of size txt, where t is the length of us
%
len = size(us,2);
Kmat = zeros(len,len);
for i = 1:len
    for j = i:len
        Kmat(i,j) = kernel(us(:,i), us(:,j));
    end
end
out = Kmat + Kmat' - diag(diag(Kmat));

end
