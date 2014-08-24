function out = ker_matrix(us,kernel)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructs the kernel matrix all at once.
% Only used to generate test functions.
% INPUTS: 
% us - control points 
% kernel - handle for the kernel generating function
% OUTPUTS
% out - matrix of size txt, where t is the length of us
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

len = size(us,2);
Kmat = zeros(len,len);
for i = 1:len
    for j = i:len
        Kmat(i,j) = kernel(us(:,i), us(:,j));
    end
end
out = Kmat + Kmat' - diag(diag(Kmat));
