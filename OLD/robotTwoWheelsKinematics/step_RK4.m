function x1 = step_RK4(h,x0,u,CON,PAR,fun)

% simulates the trajectory of all states using classical Runge-Kutta (RK4)
% the dynamics used is that given by the function handle fun

% extract constraints
x_cnstr = CON.x_cnstr;
y_cnstr = CON.y_cnstr;
w1_cnstr = CON.w1_cnstr;
w2_cnstr = CON.w2_cnstr;

% using classical Runge-Kutta method (RK4)
% assuming u is constant
k1 = h * fun(0,x0,PAR,u,false);
% do linear interpolation to get betwixt-u
x_k1 = x0 + k1/2;
k2 = h * fun(0,x_k1,PAR,u,false);
x_k2 = x0+ k2/2;
k3 = h * fun(0,x_k2,PAR,u,false);
x_k3 = x0 + k3;
k4 = h * fun(0,x_k3,PAR,u,false);
x1 = x0 + (k1 + 2*k2 + 2*k3 + k4)/6;
% contraint checking
if (abs(x1(1)) > x_cnstr || abs(x1(2)) > y_cnstr || ...
    abs(u(1)) > w1_cnstr || abs(u(2)) > w2_cnstr)
    error('Constraints violated. Aborting...');
end

