function x_dot = springDamperNominalDynamics(~,x,PAR,u,false)
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
    A = PAR.A;
    B = PAR.B;
    x_dot = A*x + B*u;
end