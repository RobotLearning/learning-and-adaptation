function out = ker_vector(t,x,us,fun)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct the kernel vector - covariances of the test pt and previous pts
% INPUTS: 
% t - last time stage
% x - test point
% us - control points to be tried (up to horizon)
% fun - function handle for the kernel
% OUTPUT
% out - kernel column vector at the previous points and current test pt.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Kvect = zeros(t,1);
for i = 1:t
    Kvect(i) = fun(us(:,i),x(:));
end
out = Kvect;
end