% see if examples have costs that match with the model

k = PAR.k;
c = PAR.c;
m = PAR.m;
dbg_cost = zeros(1,num_tot);

for t = 1:num_tot
    % assuming weight matrix is eye(2)
    x_t = CxU(t,1);
    x_tdot = CxU(t,2);
    dx_dot = CxU(t,4) - CxU(t,2);
    u = CxU(t,5);
    delta_k = (1/100) * (x_t ^3);
    dbg_cost(t) = (h^2)*(x_t^2)*(delta_k^2)/(m^2) + ...
                  (2*h/m)*x_t*delta_k*dx_dot + ...
                  (2*(h^2)/(m^2))*x_t*delta_k*(u - k*x_t - c*x_tdot);
   
end
%num_tot = 100;
%plot(1:num_tot, cost(1:num_tot), '-r', 1:num_tot, dbg_cost(1:num_tot), '*b');

plot(1:num_tot, cost, '-r', 1:num_tot, dbg_cost, '*b');