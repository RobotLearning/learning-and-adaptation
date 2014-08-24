function Kmat = ker_matrix_iter(t,us,fun,mat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct the kernel matrix iteratively
% INPUTS: 
% t - last time stage (up to horizon)
% us - control law points 
% fun - function handle for the kernel
% mat - previous kernel matrix
% OUTPUT
% Kmat - updated kernel matrix at the previous points and current test pt.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% updates iteratively
Kmat = mat; 

for i = 1:t-1
    Kmat(i,t) = fun(us(:,t), us(:,i));
    % since kernel matrix is symmetric simply update the last row with 
    % the last column
    Kmat(t,i) = Kmat(i,t);
end
% end point
Kmat(t,t) = fun(us(:,t), us(:,t)); % normally is equal to 1