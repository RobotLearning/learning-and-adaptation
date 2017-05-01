function out = kernel(x1,x2,L,type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kernel called by the ker_matrix, ker_matrix_iter and ker_vector
% functions.
% x1 - Vector 1
% x2 - Vector 2 
% L - lengthscale parameter(s) for the kernel
% type - type of the kernel: Gaussian or linear
%
% out - covariance between x1 and x2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch type
    case 'squared exponential iso'
        if length(L) == 1
            out =  exp(-(norm(x1(:)-x2(:),2)^2)/(2*(L^2)));
        else
            error('Lengthscale parameter should be scalar!');
        end
    case 'linear iso'
        if isempty(L)
            out =  x1(:)'*x2(:);
        else
            error('Lengthscale parameter should be empty!');
        end
    case 'squared exponential ard'
        InvGamma = diag(1./(L.^2));
        out = exp(-0.5*((x1(:)-x2(:))')*InvGamma*(x1(:)-x2(:)));
    case 'linear ard'
        InvGamma = diag(1./(L.^2));
        out = x1(:)'*InvGamma*x2(:);
    otherwise
        error('Unrecognized kernel type');
end

return