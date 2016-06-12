% Produce duplication and elimination matrices
% for transforming n-dim vech to vec

function [D,E] = genDuplicationMatrix(N,n,m)

% selection matrix
T = tril(ones(N));
% replace each one with nxm ones matrix 
T = kron(T,ones(n,m));
% selection vector
Tvec = T(:);
% duplication matrix
D = diag(Tvec);
D = D(logical(Tvec),:)';
% elimination matrix is its transpose
E = D';
