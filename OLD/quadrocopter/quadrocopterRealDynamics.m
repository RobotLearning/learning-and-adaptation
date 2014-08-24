function [x_dot, Aout, Bout] = quadrocopterRealDynamics(~,x,PAR,u,flg_bst)
    % [x_dot, Aout, Bout] = quadrocopterNominalDynamics(t,x,param,u,flg_bst)
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
    % used for actuator mismatch
    k = 1;
    
    %%%%%%%%%%%%%%%%% WIND DISTURBANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Disturbance caused by Wind 
    %%{
    F_y = PAR.Fan.BlowPressure * PAR.Quad.A * ...
        sin(PAR.Fan.Angle + x(5)) * cos(PAR.Fan.Angle);
    F_z = PAR.Fan.BlowPressure * PAR.Quad.A * ...
        sin(PAR.Fan.Angle + x(5)) * sin(PAR.Fan.Angle);
    F_dist = [0; F_y; 0; F_z; 0];
    %}
    
    %%%%%%%%%%%%%%%%% ACTUATOR MISMATCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    k = 1.1;
    F_dist = zeros(length(x),1);
    %}
    
    %%%%%%%%%%%%%%% DYNAMICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % x_dot = Ax + B(x)u + C
    A = [0 1 0 0 0;
         zeros(1,5);
         0 0 0 1 0;
         zeros(2,5);];
    B = [0            0;
         -k*sin(x(5)) 0;
         0            0;
         k*cos(x(5))  0;
         0            1];
    C = [0; 0; 0; -g; 0];
    
    x_dot = A*x + B*u + C + F_dist;
    
    %%%%%%%%%%%%%%% PASSING A,B MATRICES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(flg_bst)
        Aout = A*x + C;
        Bout = B;
    end
end