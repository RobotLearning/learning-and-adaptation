function x = sim_RK4(t,x_init,u_trj,~,param,fun)

% simulates the trajectory of all states using classical Runge-Kutta (RK4)
% the dynamics used is that given by the function handle fun

N = length(t) - 1;
h = t(2) - t(1);
x = zeros(length(x_init),N+1);
x(:,1) = x_init(:);

for i = 1:N
    % get trajectory of states
    % using classical Runge-Kutta method (RK4)
    k1 = h * fun(i,x(:,i),param,u_trj(:,i),false);
    x_k1 = x(:,i) + k1/2;
    u_k = interp1(t, u_trj', (t(i)+t(i+1))/2, 'linear');
    u_k = u_k';
    k2 = h * fun(i+1/2,x_k1,param,u_k,false);
    x_k2 = x(:,i) + k2/2;
    k3 = h * fun(i+1/2,x_k2,param,u_k,false);
    x_k3 = x(:,i) + k3;
    k4 = h * fun(i+1,x_k3,param,u_trj(:,i+1),false);
    x(:,i+1) = x(:,i) + (k1 + 2*k2 + 2*k3 + k4)/6;
    % no contraint checking
end