% get initial conditions from generated trajectory
x_init = x_nom(:,1);

% plot desired trajectory and deviation if nominal input is applied
figure;
title('Desired trj and deviation with nominal input');
xy_nom = x_nom(dims_tr,:);
xy_real = x_real(dims_tr,:);
for j = 1:length(dims_tr)
    subplot(2,2,j); 
    plot(t, xy_nom(j,:), '-r', t, xy_real(j,:), '-g');
    legend('nominal (desired) trj', 'trj followed with nominal u');
    hold on;
end
subplot(2,2,[3 4]); 
plot(xy_nom(1,:), xy_nom(2,:), '-r', xy_real(1,:), xy_real(2,:), '-g');
legend('nominal (desired) trj', 'trj followed with nominal u');
hold on;