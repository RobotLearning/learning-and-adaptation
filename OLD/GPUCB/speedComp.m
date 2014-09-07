%% Speed comparison between incremental inverse and backslash
clc; clear; close all;

horizon = 500;
trials  = 5;
var_noise = 0.25;
% specify a kernel
l = 0.1;
type = 'squared exponential iso';
kern = @(x1, x2) kernel(x1,x2,l,type);
u_past = rand(1,horizon);

time_incr = zeros(trials, horizon);
time_batch = zeros(trials, horizon);

for i = 1:trials
    
    Kinv = 1/(1 + var_noise);
    K = 0;
    u_past = rand(1,horizon);
    for t = 1:horizon

        x = rand;
        % calculate kernel vector k and store
        vec = ker_vector(t,x,u_past(1:t),kern);

        tic;
        if (t > 1) % if t is 1, then Kinv is a scalar
            % calculate A and store, calculate s
            Q = ker_vector(t-1,u_past(t),u_past(1:t-1),kern);
            A = Kinv * Q;
            S = kern(u_past(t), u_past(t)) + var_noise; 
            % calculate M
            m = 1/(S - Q'*A);
            % calculate Phat
            P_hat = Kinv + m*(A*A');
            % calculate Qhat
            Q_hat = -m*A;
            % form the full inverse matrix
            Kinv = [P_hat, Q_hat; Q_hat', m];
        end

        out = Kinv * vec;
        time_incr(i,t) = toc;
        
        % measure batch inverse multiplication
        tic;
        K = ker_matrix_iter(t,u_past(1:t),kern,K) + var_noise * eye(t);
        out2 = K \ vec;
        time_batch(i,t) = toc;
    end
    
    % display when slow (i.e. large horizon length)
    fprintf('Trial number %d \n', i);
end

%% Plot the running times

n = 100; % start plotting with size n
% average the results
time_incr = sum(time_incr)/trials;
time_batch = sum(time_batch)/trials;
figure(1);
plot(n:horizon, time_incr(n:horizon), '-b', ...
     n:horizon, time_batch(n:horizon), '-r');
legend('incremental inverse mult.', 'batch inverse mult.');


figure(2);
time_diff = time_batch - time_incr;
plot(n:horizon, time_diff(n:horizon), '-b');
title('\Delta Running Times');