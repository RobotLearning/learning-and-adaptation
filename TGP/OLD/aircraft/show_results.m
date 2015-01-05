mpc_horizon = STR.MPC.HORIZON;

% display SSE error
err_total_GP(run) = sum(sum((x_GP(dims_tr,:) - x_nom(dims_tr,:)).^2)); 
err_total_MPC(run) = sum(sum((x_MPC(dims_tr,:) - x_nom(dims_tr,:)).^2));
fprintf('SSE for CGP-UCB is %f \n', err_total_GP(run));
fprintf('SSE for NMPC, h = %d is %f \n', ...
        mpc_horizon, err_total_MPC(run));

% plot results
figure(2*run);
plot(x_nom(dims_tr(1),:), x_nom(dims_tr(2),:), '-r', ...
     x_MPC(dims_tr(1),:), x_MPC(dims_tr(2),:), '-k', ...
     x_GP(dims_tr(1),:), x_GP(dims_tr(2),:), '-b');
mpc_legend = sprintf('NMPC, h = %d', mpc_horizon);
legend('nominal (desired) trj', mpc_legend,'CGPUCB');

figure(2*run+1);
subplot(2,1,1);
time = h * (1:iter);
plot(time, u_trj(1,1:end-1), '-r', ...
     time, u_mpc(1,1:end-1), '-k', ...
     time, u(1,1:end-1),'-b');
legend('nominal', mpc_legend, 'CGPUCB');
subplot(2,1,2);
plot(time, u_trj(2,1:end-1), '-r', ...
     time, u_mpc(2,1:end-1), '-k', ...
     time, u(2,1:end-1),'-b');
legend('nominal', mpc_legend, 'CGPUCB');