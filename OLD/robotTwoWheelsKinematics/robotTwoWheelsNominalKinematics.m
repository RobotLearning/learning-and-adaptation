function [x_dot, Jx, Ju] = robotTwoWheelsNominalKinematics(~,x,PAR,u,flg_jcb)
    % [x_dot, Jx, Ju] = robotTwoWheelsNominalDynamics(t,x,param,u,flg_jcb)
    % differential equations of the two-wheeled robot kinematics
    % ---------------------------------------------------------------------
    % INPUTS
    %     t               time
    %     x               state vector
    %     param           model parameters
    %     u               controller input
    %     flg_jcb         flag for calculating the jacobian
    % OUTPUTS
    %     x_dot           f(x(t),u(t),t)
    %     J               Jacobian matrix
    % ---------------------------------------------------------------------
    R_1 = PAR.R1;
    R_2 = PAR.R2;
    d = PAR.d;
    % Kinematic equations
    A = zeros(length(x),length(x));
    B = [R_1/2 * cos(x(3)), R_2/2 * cos(x(3));
         R_1/2 * sin(x(3)), R_2/2 * sin(x(3));
         R_1/d,            -R_2/d];
    x_dot = A*x + B*u;
    % Jacobian calculation
    % derivatives used in the jacobian 
    if(flg_jcb)
        df1dx3 = [-R_1/2 * sin(x(3)), -R_2/2 * sin(x(3))] * u;
        df2dx3 = [R_1/2 * cos(x(3)), R_2/2 * cos(x(3))] * u;
        Jx = [0 0 df1dx3;
              0 0 df2dx3;
              0 0      0];       
        Ju = B;
    end
end