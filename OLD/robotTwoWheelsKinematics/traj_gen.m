%% Create desired nominal trajectory

% create a pre-nominal sine-curve 
x(1,:) = t;
x(2,:) = sin(2*pi*t);
x(3,1:end-1) = atan2(diff(x(2,:)),diff(x(1,:)));
x(3,end) = x(3,1); % since trajectory is periodical after t = 1
u_dim = 2;
u_trj = zeros(u_dim,N);
% maximum allowed phi difference in one discrete time step
% TODO: this is not being used!
CON.phi_dev_max = 0.05;
% initial condition for phi
phi_prev = x(3,1);

R_1 = PAR.R1;
R_2 = PAR.R2;
d = PAR.d;

for i = 1:Nu
    % discretize B(phi) around phi0
    B = h * [R_1/2 * cos(x(3,i)), R_2/2 * cos(x(3,i));
             R_1/2 * sin(x(3,i)), R_2/2 * sin(x(3,i));
             R_1/d,                 -R_2/d];
    % deviation along desired trajectory
    delta = x(2:3,i+1) - x(2:3,i);
    % solve for u
    u_trj(:,i) = B(2:3,:)\delta;
end

u_trj(:,end) = u_trj(:,end-1);

% simulate nominal trajectory with RK4
handle = @robotTwoWheelsNominalKinematics;
x = sim_RK4(t,x,u_trj,CON,PAR,handle);

% interpolate to plot till N
% plot trajectories
figure(1);
plot(t, x);
legend('robot displacement x', 'robot displacement y', ...
    'angular yaw phi');
title('Nominal trajectory to be followed');
figure(2);
subplot(2,3,1);
plot(t, x(1,:),'r'); hold on;
title('x of the trajectory over iterations');
subplot(2,3,2);
plot(t, x(2,:),'r'); hold on;
title('y of the trajectory over iterations');
subplot(2,3,3);
plot(t, x(3,:),'r'); hold on;
title('phi of the trajectory over iterations');
subplot(2,3,4);
plot(t, u_trj(1,:),'r'); hold on;
title('Control input w1 over iterations');
subplot(2,3,5);
plot(t, u_trj(2,:),'r'); hold on;
title('Control input w2 over iterations');