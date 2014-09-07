function [t,u_trj,x] = traj_gen(trj,PAR,h)

g = PAR.g;
m = PAR.m;
J = PAR.J;
r = PAR.r;
N = length(trj);
t = linspace(h,h*N,N);

z1 = trj(1,:);
z2 = trj(2,:);
ddz1 = trj(3,:);
ddz2 = trj(4,:);

theta = atan(-ddz1 ./ (ddz2 + g));
xt = z1 + (J/(m*r))*sin(theta);
yt = z2 - (J/(m*r))*cos(theta);
dxt = diff(xt)./diff(t);
ddxt = diff(dxt)./diff(t(1:end-1));
dyt = diff(yt)./diff(t);
dtheta = diff(theta)./diff(t);
ddtheta = diff(dtheta)./diff(t(1:end-1));

u_trj(1,:) = ddtheta * J/r;
u_trj(2,:) = (u_trj(1,:).*cos(theta(1:end-2)) - (m * ddxt))./ ...
              sin(theta(1:end-2)); 

u_trj(:,N-1) = u_trj(:,N-2);
u_trj(:,N) = u_trj(:,N-1);

% simulate nominal trajectory with RK4
handle = @aircraftNominalDynamics;
x_init = [xt(1); dxt(1); yt(1); dyt(1); theta(1); dtheta(1)];
x = sim_RK4(t,x_init,u_trj,0,PAR,handle);

end