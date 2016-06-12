%% Estimates the whole LTV block matrices together

% inputs u and e are the inputs and errors (outputs)
% where
%
% u is mN x K,
% e is nN x K,
%
% K = # of experiments
% N = horizon size
% m = number of inputs
% n = number of outputs
% 
% Estimating the matrices in the form [B1; A2*B1; A3*A2*B1; ...]
% and then transforming back to block-Hankel form
%
function F_est = estimateSystemMatrices(u,e,systemSize)

    n = systemSize(1);
    m = systemSize(2);
    K = size(u,2);
    assert(size(u,1)/m == size(e,1)/m, 'size of horizon does not match!');
    N = size(u,1)/m;
    F_est = zeros(n*N,m*N);

    % estimate F
    for i = 1:K
        U((N*(i-1))+1:(N*(i-1))+N,1:m*N*N) = kron(eye(N),u(:,i)');
        E((N*(i-1))+1:(N*(i-1))+N,1:n) = reshape(e(:,i)',n,N)';
    end
    
    D = duplicateVech(N,m);
    M = U * D';
    Est = pinv(M) * E;
    EstAddZero = Est' * D;

    for i = 1:N
        F_est((i-1)*n+1:i*n,:) = EstAddZero(:,N*m*(i-1)+1:N*m*i);
    end    


end