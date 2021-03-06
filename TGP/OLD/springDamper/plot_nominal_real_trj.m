% get initial conditions from generated trajectory
x_init = x_nom(:,1);

% plot desired trajectory and deviation if nominal input is applied
figure;
title('Desired trj and deviation with nominal input');
yz_nom = x_nom(dims_tr,:);
yz_real = x_real(dims_tr,:);
for j = 1:dims_tr
    subplot(2,1,j); 
    plot(t, yz_nom(j,:), '-r', t, yz_real(j,:), '-g');
    legend('nominal (desired) trj', 'trj followed with nominal u');
    hold on;
end
subplot(2,1,dims_tr+1); 
plot(t, u_trj, '-r');
legend('nominal u');
hold on;