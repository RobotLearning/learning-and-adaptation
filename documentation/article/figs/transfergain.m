clear; clc;
ilc = [0.534257, 2.211294, 2.305692];
tgp = [2.552392, 3.562660, 3.837903];
tgps = [0.222369, 1.007975, 0.543684];

traj = 1:length(ilc);

plot(traj,ilc,'x', traj,tgp, '*', traj, tgps, 'o');
axis([0.5 length(ilc)+0.5  0 4]);
set(gca,'XTick',1:length(ilc));
l = legend('ILC','TGP','TGP-with prior' );
set(l,'Location','NorthWest');
xlabel('Trajectories');
ylabel('Sum of Squares Error (SSE)');