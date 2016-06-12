% Function to test ILC methods here
% The idea is to check which one is better in terms of 
% Monotonicity when model mismatch is present

function u = ilc(U,E,alg,model)

    F = model;
    switch alg
        % NO MODEL LEARNING
        % FIRST ORDER ILC
        case 'simple'            
            assert(size(F,1) == size(F,2),...
                'Cannot apply PD-ILC, matrix is not square!');
            alpha = 0.2;
            u_ilc = alpha * E(:,end);
            u = U(:,end) + u_ilc;
        % THESE ARE ALL MODEL BASED
        % SECOND ORDER - PLANT INVERSION - BASED METHODS
        case 'pinv'
            u_ilc = pinv(F) * E(:,end);
            u = U(:,end) - u_ilc;
        % PSEUDOINVERSE TRUNCATION CORRESPONDS TO PCA (95% limit)
        case 'pinv-truncate'
            u_ilc = pinv(F,0.05) * E(:,end);
            u = U(:,end) - u_ilc;
        case 'inv'
            u_ilc = F \ E(:,end);
            u = U(:,end) - u_ilc;
        % REGULARIZATION CORRESPONDS TO PENALIZING CONTROL INPUT CHANGEs
        case 'regular-inv'
            R = 0.4 * eye(size(F,2));
            u_ilc = (F'*F + R)\(F'*E(:,end));
            u = U(:,end) - u_ilc;
            
        % FIRST ORDER
        % MORE RESISTANT TO MODEL MISSPECIFICATION
        case 'grad-descent'
            u_ilc = F'*E(:,end);
            u = U(:,end) - u_ilc;
            
        otherwise
            error('ILC type not recognized!');
    end

end