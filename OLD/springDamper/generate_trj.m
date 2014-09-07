% determine coordinates of the new wave
A_d = 10;
w_d = 4;
trj = A_d * sin(w_d * time);
dtrj = w_d*2*pi * A_d * cos(w_d*2*pi*time);
traj = [trj; dtrj];

[time,u_trj,x_nom] = traj_gen(traj,PAR,h);
N = length(time);