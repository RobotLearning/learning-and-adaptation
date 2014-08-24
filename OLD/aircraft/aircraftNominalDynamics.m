function x_dot = aircraftNominalDynamics(~,x,PAR,u,~)
    % x_dot = aircraftrNominalDynamics(t,x,param,u,flg_jcb)
    % differential equations of the quadrocopter
    % ---------------------------------------------------------------------
    % INPUTS
    %     t               time
    %     x               state vector
    %     param           model parameters
    %     u               controller input
    % OUTPUTS
    %     x_dot           f(x(t),u(t),t)
    % ---------------------------------------------------------------------
    
    % Extract parameters
    g = PAR.g; 
    r = PAR.r;
    J = PAR.J;
    m = PAR.m;
    % Dynamics
    % x_dot = Ax + B(x)u + C
    A = [0 1 0 0 0 0;
         zeros(1,6);
         0 0 0 1 0 0;
         zeros(1,6);
         0 0 0 0 0 1;
         zeros(1,6)];
    B = [0            0;
         cos(x(5))/m  -sin(x(5))/m;
         0            0;
         sin(x(5))/m  cos(x(5))/m;
         0            0;
         r/J          0];
    C = [0; 0; 0; -g; 0; 0];
    x_dot = A*x + B*u + C;
end