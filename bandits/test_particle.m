%% Testing simple gradient based learning

% Testing cautious learning using gradients
clc; clear; close all;

%% Plot a function

rng(2);
hp.type = 'squared exponential ard';
hp.l = 0.25;
hp.scale = 1;
hp.noise.var = 0.0001;

n = 50;
mu = 0;
x = linspace(0,1,n)';
l = hp.l;
s = hp.scale;
kern = @(x1,x2) s * exp(-(norm(x1(:)-x2(:),2)^2)/(2*(l^2)));
Sigma = zeros(n,n);
for i = 1:n
    for j = 1:n
        Sigma(i,j) = kern(x(i),x(j));
    end
end
[U,S] = eig(Sigma);
f = U * sqrt(max(0,S)) * randn(n,1);

plot(x,f);
hold on;

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

% Fit a GP with 5 noisy samples
n_sample = 5;
idx_sample = randi(n,n_sample,1);
x_sample = x(idx_sample);
y_sample = f(idx_sample) + sqrt(hp.noise.var)*randn(n_sample,1);
gp = GP(hp,x_sample',y_sample,true);

m = zeros(n,1);
s2 = zeros(n,1);
m_der = zeros(n,1);
s2_der = zeros(n,1);
for i = 1:n
    [mean, var, mean_der, var_der] = gp.predict(x(i));
    m(i) = mean;
    s2(i) = var;
    m_der(i) = mean_der;
    s2_der(i) = var_der;
end
plot(x,m,'r');
f = [m+2*sqrt(s2); flip(m-2*sqrt(s2))];
fill([x;flip(x)],f,'b','FaceAlpha',0.3);
plot(x_sample, y_sample, '*', 'Color', 'k', 'MarkerSize', 8);

%% Climb up to local max of function with a particle

% solve bvp starting from x0 = 0 and ends up at time T to f_max
%[f_max,idx] = max(f_poly);
[x_ext_cand,~] = fsolve(@(x)gp.mean_der(x),0);
xf = x_ext_cand;
x0 = 0.0;
T = 5.0;
fnc_particle = @(t,x) [x(2); -gp.mean_der(x(1))];
f_bc = @(x1,x2) [x1(1); x2(1) - xf];
n_init = 10;
x_init = @(t) [x0 + (t/T) * (xf-x0); 0.1];
solinit = bvpinit(linspace(0,T,n_init),x_init);
sol = bvp4c(fnc_particle,f_bc,solinit);
x_particle = sol.y(1,:);
f_particle = gp.predict_mesh(x_particle);
scatter(x_particle,f_particle);
% 
% %% Climb down to go to another max
% 
% x_last = sol.y(:,end);
% [t,x_cnt] = ode45(fnc_particle,[0 5],x_last);
% f_cnt = f_model(x_cnt(:,1));
% scatter(x_cnt(:,1),f_cnt);