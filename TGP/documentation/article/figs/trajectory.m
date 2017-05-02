figure
subplot(1,3,1);
plot(trj.s(3,:),trj.s(1,:),'--','color','black'); hold on;
grayShades = [0.8,0.7,0.6,0.4,0.2,0.1];
for i = 1:5
    plot(trj.PERF(1+i).x(3,:),trj.PERF(1+i).x(1,:),'color', ones(1,3)*grayShades(i))
end 
plot(trj.PERF(13).x(3,:),trj.PERF(13).x(1,:),'color', ones(1,3)*grayShades(6))
axis([0 0.3 0 0.35]);
hold off;
title('TGP');

subplot(1,3,2);
plot(trj.s(3,:),trj.s(1,:),'--','color','black'); hold on;

for i = 1:5
    plot(trj.PERF(6+i).x(3,:),trj.PERF(6+i).x(1,:),'color', ones(1,3)*grayShades(i))
end 
plot(trj.PERF(14).x(14),trj.PERF(6+i).x(14),'color', ones(1,3)*grayShades(6))
axis([0 0.3 0 0.35]);
hold off;
title('ILC');

subplot(1,3,3);
plot(trj.s(3,:),trj.s(1,:),'--','color','black'); hold on;
for i = 1:1
    plot(trj.PERF(11+i).x(3,:),trj.PERF(11+i).x(1,:),'color', ones(1,3)*grayShades(i))
end 
axis([0 0.3 0 0.35]);
hold off;
title('MPC')