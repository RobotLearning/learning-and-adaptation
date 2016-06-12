%% Recursively estimate LTV block system matrices using LBR

% inputs u and e are the inputs and errors (outputs)
% where
%
% u is mN x 1,
% e is nN x 1,
%
% N = horizon size
% m = number of inputs
% n = number of outputs
% 
% Estimating the matrices in vectorized form
%
function [F_est,Gamma,hp] = recursiveEstSysMat(F_est,Gamma,hp,u,e,N)

    m = length(u)/N;
    n = length(e)/N;
    [Est,D,~] = strip(F_est,m,n,N);

    % estimate F
    M = kron(u',eye(n*N))*D;
    %Est = pinv(M)*e(:);
    %EstAddZero = D*Est;
    
    [Est,Gamma,hp] = LBR(Est,Gamma,M,e,hp);
    EstAddZero = D*Est;
    F_est = reshape(EstAddZero',n*N,m*N);

end

% Strip the previously estimated model matrix
function [Est,D,E] = strip(F_est,m,n,N)

    Vec = F_est(:);
    [D,E] = genDuplicationMatrix(N,n,m);
    Est = E*Vec;
    
end