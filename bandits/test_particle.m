%% Testing simple gradient based learning

% Testing cautious learning using gradients
clc; clear; close all;

%% Plot a function

rng(3);
hp.type = 'squared exponential ard';
hp.l = 0.2;
hp.scale = 1;
hp.noise.var = 0.0001;

n = 50;
mu = 0;
mesh = linspace(0,1,n)';
l = hp.l;
s = hp.scale;
kern = @(x1,x2) s * exp(-(norm(x1(:)-x2(:),2)^2)/(2*(l^2)));
Sigma = ker_matrix(mesh',kern);
[U,S] = eig(Sigma);
f = U * sqrt(max(0,S)) * randn(n,1);

plot(mesh,f);
hold on;
legend('function');

%% Traverse the path with Thompson sampling
%{
horizon = 20;
buffer = 20;
cand_max = zeros(buffer,1);
cand_min = zeros(buffer,1);
gp = GP(hp,[],[]);
x = mesh(20);
plot(x, f(20), '*', 'Color', 'k', 'MarkerSize', 8);
% start from 0
for i = 1:horizon
    pause(0.5);
    %input('Press return to continue\n');
    [mu,Sigma] = gp.predict_mesh(mesh');
    [U,S] = eig(Sigma);
    for j = 1:buffer
        f_thompson = mu(:) + U * sqrt(max(0,S)) * randn(n,1);
        [~,cand_max(j)] = max(f_thompson);
        [~,cand_min(j)] = min(f_thompson);
    end
    cand_min(buffer+1:buffer+2) = [1,n];
    
    % calculate force
    force = sum(1/(mesh(cand_max) - x)) - ...
        sum(1/(mesh(cand_min) - x));
    x = min(1,max(0,x + force/10));
    %plot(mesh,f_thompson);
    y_sample = interp1(mesh,f,x) + sqrt(hp.noise.var)*randn;
    gp.update(x,y_sample);
    plot(x, y_sample, '*', 'Color', 'k', 'MarkerSize', 8);
end
hold off;
%}
    
%% OTHER CODE
% Fit a fourth order model
%{
phi = @(x) [ones(length(x),1), x, x.^2, x.^3, x.^4];
phi_deriv = @(x) [zeros(length(x),1), 1, 2*x, 3*x.^2, 4*x.^3];
Phi = phi(x);
a = Phi \ f;
f_poly = Phi*a;
f_model = @(x) phi(x)*a;
f_model_deriv = @(x) phi_deriv(x)*a;
[x_ext_cand,~] = fsolve(f_model_deriv,0);
hold on;
plot(x,f_poly);
legend('function', 'model surrogate');
%}

% Gradient Descent with EI utility function
%{
% n_sample = 2;
% idx_sample = 1:n_sample; %randi(n,n_sample,1);
% x_sample = x(idx_sample);
% y_sample = f(idx_sample) + sqrt(hp.noise.var)*randn(n_sample,1);
% gp = GP(hp,x_sample',y_sample,true);
% 
% m = zeros(n,1);
% s2 = zeros(n,1);
% m_der = zeros(n,1);
% s2_der = zeros(n,1);
% [m,s2] = gp.predict_mesh(x');
% s2 = diag(s2);
% plot(x,m,'r');
% f_ucb = [m+2*sqrt(s2); flip(m-2*sqrt(s2))];
% fill([x;flip(x)],f_ucb,'b','FaceAlpha',0.3);
% plot(x_sample, y_sample, '*', 'Color', 'k', 'MarkerSize', 8);
% 
% density = @(x) (1/sqrt(2*pi)) .* exp((-1/2).*(x.^2));
% distr = @(x) (1 - erf(x./sqrt(2)))./2;
% acq = @(m,s2) (m - max(m)) .* distr((max(m) - m)./sqrt(s2)) + ...
%           sqrt(s2) .* density((max(m) - m)./sqrt(s2));
% [mu,Sigma] = gp.predict_mesh(x');
% s2 = diag(Sigma);
% ei = 10 * acq(mu,s2);
% plot(x,ei,'--k');
% legend('fnc','mean','mean+2*var','samples','ei');
%}

% Climb up to local max of function with a particle
%{
% solve bvp starting from x0 = 0 and ends up at time T to f_max
%[f_max,idx] = max(f_poly);
% acq_func = @(x) gp.mean(x) + 0.5 * gp.var(x);
% acq_der = @(x) gp.mean_der(x) + 0.5 * gp.var_der(x);
% %[x_ext_cand,~] = fsolve(acq_der,0.5);
% [~,idx_cand] = max(acq_func(x));
% xf = x_ext_cand;
% x0 = 0.0;
% T = 5.0;
% fnc_particle = @(t,x) [x(2); -acq_der(x(1))];
% f_bc = @(x1,x2) [x1(1); x2(1) - xf];
% n_init = 10;
% x_init = @(t) [x0 + (t/T) * (xf-x0); 0.1];
% solinit = bvpinit(linspace(0,T,n_init),x_init);
% sol = bvp4c(fnc_particle,f_bc,solinit);
% x_particle = sol.y(1,:);
% f_particle = interp1(x,f,x_particle);
% scatter(x_particle,f_particle);
% 
% %% Climb down to go to another max
% 
% x_last = sol.y(:,end);
% [t,x_cnt] = ode45(fnc_particle,[0 5],x_last);
% f_cnt = f_model(x_cnt(:,1));
% scatter(x_cnt(:,1),f_cnt);
%}