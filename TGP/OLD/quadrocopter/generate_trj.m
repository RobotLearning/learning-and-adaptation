traj = 'Wave';
% determine coordinates of the new wave
%y_coord = [0 0.5 2.5 3.0]/10;
%z_coord = [0 1.2 0.0 1.2]/10;

incr_y = [0.126, 0.145, 0.309];
incr_z = [0.1833, 0.412, 0.4691];
y_coord = ([0 1.2 2.4 3.2] + [0 incr_y])/10;
z_coord = ([0 1.8 0.0 1.2] + [0 incr_z])/10;

%{
alpha = pi/2;
rot = [cos(alpha), sin(alpha);
       -sin(alpha), cos(alpha)];
newWave = rot * [y_coord; z_coord];
y_coord = newWave(1,:);
z_coord = newWave(2,:);
%}

[time,u_trj,x_nom] = quad_traj_gen(traj,PAR,CON,h,y_coord,z_coord,1);
N = length(time);