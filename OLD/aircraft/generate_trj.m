% determine coordinates of the new wave
% constructing z1 and z2
r = 10;
w = pi/5;
t = linspace(h,h*N,N);

z1 = r - r*cos(w*t); 
z2 = r*sin(w*t);
ddz1 = (w^2)*r*cos(w*t);
ddz2 = -(w^2)*r*sin(w*t);
traj = [z1; z2; ddz1; ddz2];
[t,u_trj,x_nom] = traj_gen(traj,PAR,h);

for j = 1:length(dims_tr)
    subplot(2,2,j); 
    plot(t, x_nom(j,:), '-r');
    legend('nominal (desired) trj');
    %hold on;
end

subplot(2,2,[3 4]); 
plot(x_nom(1,:), x_nom(2,:), '-r');
legend('nominal (desired) trj');
%hold on;