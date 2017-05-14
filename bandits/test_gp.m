%% Test Gaussian Process (GP) regression

clc; clear; close all;

%% Simple regression
%rng('default');
hp.type = 'squared exponential ard';
hp.l = 0.25;
hp.scale = 1;
hp.noise.var = 0.01;

n = 20;
mu = 0;
x = randn(n,1);
l = hp.l;
s = hp.scale;
kern = @(x1,x2) s * exp(-(norm(x1(:)-x2(:),2)^2)/(2*(l^2)));
Sigma = zeros(n,n);
for i = 1:n
    for j = 1:n
        Sigma(i,j) = kern(x(i),x(j));
    end
end
y = mu + chol(Sigma + hp.noise.var*eye(n)) * x;
gp = GP(hp,x',y);

meshsize = 100;
z = linspace(min(x) - 0.5, max(x) + 0.5, meshsize);
[m,s2] = gp.predict_mesh(z);
s2 = diag(s2);
% m = zeros(meshsize,1);
% s2 = zeros(meshsize,1);
% for i = 1:length(z)
%     [mean, var] = gp.predict(z(i));
%     m(i) = mean;
%     s2(i) = var;
% end

figure(1);
subplot(2,1,1);
title('GP Regression');
%set(gca, 'FontSize', 24);
f = [m+2*sqrt(s2); flip(m-2*sqrt(s2))];
fill([z';flip(z')],f,'b','FaceAlpha',0.3);
hold on; 
plot(z, m, 'LineWidth', 2, 'Color', 'r'); 
plot(x, y, '*', 'Color', 'k', 'MarkerSize', 8);
axis tight;
%axis([-1.9 1.9 -0.9 3.9])
%grid on
xlabel('input x');
ylabel('output y');
hold off;

%% Test derivatives
der = true;
gp = GP(hp,x',y,der);
m_der = zeros(meshsize,1);
s2_der = zeros(meshsize,1);
for i = 1:length(z)
    [mean, var, mean_der, var_der] = gp.predict(z(i));
    m_der(i) = mean_der;
    s2_der(i) = var_der;
end
subplot(2,1,2);
%title('GP derivatives');
plot(z,m_der,z,s2_der, 'LineWidth', 2);
m_diff = diff(m)./diff(z)'; m_diff(end+1) = m_diff(end);
hold on;
s2_diff = diff(s2)./diff(z)'; s2_diff(end+1) = s2_diff(end);
plot(z,m_diff,z,s2_diff, 'LineWidth', 2);
legend('mean der','var der','mean diff', 'var diff');
axis tight;
hold off;

%% Hyperparameter estimation

% gp.fit_hp(hp);
% z = linspace(-2, 2, 100);
% for i = 1:length(z)
%     [mean2, var2] = gp.predict(z(i));
%     m(i) = mean2;
%     s2(i) = var2;
% end
% 
% figure(2)
% title('GP regression + HP estimation');
% set(gca, 'FontSize', 24)
% f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)];
% fill([z; flip(z,1)], f, [7 7 7]/8);
% hold on; 
% plot(z, m, 'LineWidth', 2); 
% plot(x, y, '+', 'MarkerSize', 12)
% %axis([-1.9 1.9 -0.9 3.9])
% grid on
% xlabel('input, x')
% ylabel('output, y')
% hold off;