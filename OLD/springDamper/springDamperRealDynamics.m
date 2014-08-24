function x_dot = springDamperRealDynamics(~,x,PAR,u,~)
    % [x_dot,u] = systemDynamics(~,x,param,g,K)
    % differential equations of the nominal system dynamics
    % ---------------------------------------------------------------------
    % INPUTS
    %     t               time
    %     x               state vector
    %     u               control input
    % OUTPUTS
    %     x_dot           time derivatie of the state vector
    % ---------------------------------------------------------------------
    
    m = PAR.m; % spring mass
    c = PAR.c; % spring damper constant
    d = x(1); % displacement
    k_real = (d.^3 + 100*d)/100; % spring nonlinear param.
    %k_real = 1.1;
    A = [0, 1; -k_real/m, -c/m]; 
    % A changes to accommodate nonlinear spring
    B = PAR.B; 
    % B does not change
    
    x_dot = A*x + B*u;
end