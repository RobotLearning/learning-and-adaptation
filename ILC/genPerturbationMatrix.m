% Create a block lower-diagonal perturbation matrix 
% that has 2-norm (max singular value) less than a given positive alpha

function E = genPerturbationMatrix(n,m,N,a)

    assert(a > 0, 'alpha should be positive!');

    mat = randn(n*N,m*N);
    matlow = tril(mat);
    s = max(svd(matlow));
    scale = ceil(s/a);
    E = matlow/scale;
    
end