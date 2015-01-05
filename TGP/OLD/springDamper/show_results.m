mpc_horizon = STR.MPC.HORIZON;

% display SSE error
err_total_GP(run) = sum(sum((x_GP(dims_tr,:) - x_nom(dims_tr,:)).^2)); 
err_total_MPC(run) = sum(sum((x_MPC(dims_tr,:) - x_nom(dims_tr,:)).^2));
fprintf('SSE for CGP-UCB is %f \n', err_total_GP(run));
fprintf('SSE for NMPC, h = %d is %f \n', ...
        mpc_horizon, err_total_MPC(run));

% plot results
% dims_tr is a scalar in this case
figure(2*run-1);
plot(time, x_nom(dims_tr,:), '-r', ...
     time, x_MPC(dims_tr,:), '-k', ...
     time, x_GP(dims_tr,:), '-b');
mpc_legend = sprintf('NMPC, h = %d', mpc_horizon);
legend('nominal (desired) trj', mpc_legend,'CGPUCB');
title('Trajectory tracking');

figure(2*run);
subplot(2,1,1);
plot(time(1:end-1), u_trj(1,1:end-1), '-k', ...
     time(1:end-1), u_mpc(1,1:end-1),'-b');
legend('nominal', mpc_legend);
subplot(2,1,2);
plot(time(1:end-1), u_trj(1,1:end-1), '-k', ...
     time(1:end-1), u(1,1:end-1),'-b');
legend('nominal', 'CGPUCB');
title('Control inputs');