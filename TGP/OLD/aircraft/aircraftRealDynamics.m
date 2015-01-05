function x_dot = aircraftRealDynamics(~,x,PAR,u,~)
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
    
    %%%%%%%%%%%%%%%% WIND DISTURBANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Disturbance caused by Fan 
    % Define the trapezoid region with alpha and beta parameters
    %{
    F_x = PAR.Wind.Pressure * PAR.Craft.A * ...
        sin(PAR.Wind.Angle + x(5)) * cos(PAR.Wind.Angle);
    F_y = PAR.Wind.Pressure * PAR.Craft.A * ...
        sin(PAR.Wind.Angle + x(5)) * sin(PAR.Wind.Angle);
    F_dist = [0; F_x; 0; F_y; 0; 0];
    %}

    %%%%%%%%%%%%%%% GRAVITY ON PLANET X %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%{
    g_mars = 3.75;
    g_uranus = 10.5;
    g = g_uranus;
    F_dist = zeros(6,1);
    %}
    
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
    x_dot = A*x + B*u + C + F_dist/m;
end