function x1 = step_RK4(h,x0,u,~,param,fun)

% simulates the trajectory of all states using classical Runge-Kutta (RK4)
% the dynamics used is that given by the function handle fun

% using classical Runge-Kutta method (RK4)
% assuming u is constant
k1 = h * fun(0,x0,param,u,false);
x_k1 = x0 + k1/2;
k2 = h * fun(0,x_k1,param,u,false);
x_k2 = x0+ k2/2;
k3 = h * fun(0,x_k2,param,u,false);
x_k3 = x0 + k3;
k4 = h * fun(0,x_k3,param,u,false);
x1 = x0 + (k1 + 2*k2 + 2*k3 + k4)/6;
% no contraint checking
end

