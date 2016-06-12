% Linear Bayesian Regression
% updates mean and precision matrix instead of covariance
function [mu,Gamma,hp] = LBR(mu,Gamma,X,y,hp)

    n = length(mu);
    % hyperparameters of the precision matrix prior
    hp.a = hp.a + n/2;
    hp.b = hp.b + 0.5 * (y'*y + mu'*Gamma*mu);
    % Pseudoinverse seems to be more stable when there is no noise
    mu = pinv(X'*X + Gamma)*(Gamma*mu + X'*y);
    Gamma = (X'*X + Gamma);
    hp.b = hp.b - 0.5 * (mu'*Gamma*mu);
    
end