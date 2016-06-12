%% Class for testing ILC algorithms

classdef ILC < handle
    
    properties
        
        % For PD type update
        alpha; 
        % model to be considered for inversion
        model
        % parameters for LBR as a structure
        LBR
        
    end
    
    methods
        
        % Constructor for estimator
        function obj = ILC(F)
            
            obj.model = F;
            obj.alpha = 0.2;
        end
        
        % Simple PD type update rule
        function u = simplePDrule(obj,U,E)
            assert(size(E,1) == size(E,1),...
                'Cannot apply PD-ILC, matrix is not square!');
            u_ilc = obj.alpha * E(:,end);
            u = U(:,end) + u_ilc;
        end
        
        %% Model based update rules, no learning
        % Plant inversion, no learning
        % R is regularization
        % tol is truncation parameter
        function u = invertModel(obj,U,E,R,tol)
            F = obj.model;
            if R == 0
                u_ilc = pinv(F,tol) * E(:,end);
            else
                assert(size(R,1) == size(R,2),'R should be square!');
                assert(size(R,1) == size(F,2),'R should have same column length');
                u_ilc = (F'*F + R)\(F'*E(:,end));
            end
            u = U(:,end) - u_ilc;
        end
        
        % Gradient based update rule
        % More resistant to model misspecification than 2nd order methods
        function u = gradientDescent(obj,U,E)
            F = obj.model;
            u = U(:,end) - F'*E(:,end);
        end
        
        %% Adaptive ILC, learning the dynamics model
        
        % Linear Bayesian Regression
        % to estimate the mean model matrix and
        % additionally its variance
        function u = LBR(obj,U,E)
            % TODO:
        end
        
    end
    
end