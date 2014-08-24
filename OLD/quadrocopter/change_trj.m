traj = 'Wave';

% determine randomized coordinates of the new wave
%incr_y = 0.5*rand(1,3);
%incr_z = 0.5*rand(1,3);
%incr_y = [0.126, 0.145, 0.309];
%incr_z = [0.1833, 0.412, 0.4691];

%{
if run == runs - 1
y_coord = ([0 1.2 2.4 3.2] + [0 incr_y])/10;
z_coord = ([0 1.8 0.0 1.2] + [0 incr_z])/10;
else
y_coord = [0 0.5 2.5 3.0]/10;
z_coord = [0 1.2 0.0 1.2]/10;
end
%}


%%{
%rotate wave by alpha degrees
alpha = 0;

rot = [cos(alpha), sin(alpha);
       -sin(alpha), cos(alpha)];
newWave = rot * [y_coord; z_coord];
y_coord = newWave(1,:);
z_coord = newWave(2,:);
%}

[time,u_trj,x_nom] = traj_gen(traj,PAR,CON,h,y_coord,z_coord,run+1);
N = length(time);
iter = N-1;

% initialize past control signals
STR.CGP.u_past = zeros(dim,1); 
STR.MPC.u_past = zeros(dim,1);