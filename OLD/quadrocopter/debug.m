% see if examples have costs that match with the model

%% Debugging wind disturbance case

dbg_y = zeros(1,num_tot);
dbg_cost = zeros(1,num_tot);
alpha = PAR.Fan.BlowPressure * PAR.Quad.A;
q1 = Q(1,1);
q2 = Q(2,2);

%% Debugging cost differences

for t = 1:num_tot
    % assuming weight matrix is eye(6)
    phi = CxU(t,5);
    dy_dot = CxU(t,2) - CxU(t,7);
    u1 = CxU(t,11);
    u2 = CxU(t,12);
    dbg_cost(t) = h*h*(alpha*alpha - 2*alpha*u1)*(sin(phi+(h/2)*u2)^2) + ...
                  2*h*alpha*dy_dot*sin(phi+(h/2)*u2);
   
end

plot(1:num_tot, cost(1:num_tot), '-r', 1:num_tot, dbg_cost, '*b');