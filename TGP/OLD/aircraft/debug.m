% see if examples have costs that match with the model

%% Debugging gravity mismatch case
g = PAR.g;
g_uranus = 10.5;
m = PAR.m;
dbg_cost = zeros(1,num_tot);
delG = g_uranus - g;

for t = 1:num_tot
    % assuming weight matrix is eye(6)
    theta = CxU(t,5);
    dy_dot = CxU(t,10) - CxU(t,4);
    u1 = CxU(t,13);
    u2 = CxU(t,14);
    dbg_cost(t) = h*delG*(h*delG + dy_dot + (h/m)*u1*sin(theta) + ...
                          (h/m)*u2*cos(theta) - h*g);
   
end

plot(1:num_tot, cost, '-r', 1:num_tot, dbg_cost, '*b');