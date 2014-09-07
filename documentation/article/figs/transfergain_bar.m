clear; clc;
ilc = [0.534257, 2.211294, 2.305692];
tgp = [2.552392, 3.562660, 3.837903];
tgps = [0.222369, 1.007975, 0.543684];

D = [ilc',tgp',tgps'];

traj = 1:length(ilc);

hArray = bar(traj,D);%,'histic');
set(hArray(1),'EdgeColor','black','FaceColor','white');
set(hArray(2),'EdgeColor','black','FaceColor',[192,192,192]/256.);
set(hArray(3),'EdgeColor','black','FaceColor',[128,128,128]/256.);
%axis([0.5 length(ilc)+0.5  0 4]);
%set(gca,'XTick',1:length(ilc));
l = legend('ILC','TGP','TGP-with transfer' );
set(l,'Location','NorthWest');
xlabel('Trajectories');
ylabel('Sum of Squares Error (SSE)');