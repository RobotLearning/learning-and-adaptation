%% Estimates the whole LTV block matrices together (batch)

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
function F_est = batchEstSysMat(u,e,systemSize)

    n = systemSize(1);
    m = systemSize(2);
    K = size(u,2);
    assert(size(u,1)/m == size(e,1)/n, 'size of horizon does not match!');
    N = size(u,1)/m;
    F_est = zeros(n*N,m*N);

    % estimate F
    for i = 1:K
        U((N*(i-1))+1:(N*(i-1))+N,1:m*N*N) = kron(eye(N),u(:,i)');
        E((N*(i-1))+1:(N*(i-1))+N,1:n) = reshape(e(:,i)',n,N)';
    end
    
    D = duplicateBlockVech(N,m);
    M = U * D;
    Est = pinv(M) * E;
    EstAddZero = Est' * D';

    for i = 1:N
        F_est((i-1)*n+1:i*n,:) = EstAddZero(:,N*m*(i-1)+1:N*m*i);
    end    


end

% Produce duplication matrix for transforming n-dim block vech to vec

function D = duplicateBlockVech(N,m)

    % assuming A matrix is nxm

    pos = ones(1,N^2*m);
    vec = [];
    for i = 1:N-1
        val = ((i-1)*m*(N+1))+m+1:i*N*m; %i*N*n*m+1:i*(N+1)*n*m;
        vec = [vec,val]; 
        pos(val) = 0;
    end

    D = diag(pos);
    D = D(setdiff(1:N^2*m,vec),:);
    D = D';

end