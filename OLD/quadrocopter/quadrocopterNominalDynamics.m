function [x_dot, Aout, Bout] = quadrocopterNominalDynamics(~,x,PAR,u,flg_bst)
    % [x_dot, Jx, Ju] = quadrocopterNominalDynamics(t,x,param,u,flg_bst)
    % differential equations of the quadrocopter
    % ---------------------------------------------------------------------
    % INPUTS
    %     t               time
    %     x               state vector
    %     param           model parameters
    %     u               controller input
    %     flg_bst         flag for returning A,B matrices
    % OUTPUTS
    %     x_dot           f(x(t),u(t),t)
    %     Aout, Bout      returned matrices
    % ---------------------------------------------------------------------
    
    % Extract parameters
    g = PAR.g; 
    % Dynamics
    % x_dot = Ax + B(x)u + C
    A = [0 1 0 0 0;
         zeros(1,5);
         0 0 0 1 0;
         zeros(2,5);];
    B = [0          0;
         -sin(x(5)) 0;
         0          0;
         cos(x(5))  0;
         0          1];
    C = [0; 0; 0; -g; 0];
    x_dot = A*x + B*u + C;

    if(flg_bst)
        Aout = A*x + C;
        Bout = B;
    end
end