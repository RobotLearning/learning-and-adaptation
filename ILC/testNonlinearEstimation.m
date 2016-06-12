%% Testing monotonic learning control 

% 1d nonlinear system

clc; clear; close all;

eps = 0.2;
dyn = @(x,u) x^2 + eps*x - u;
nomDyn = @(x,u) x^2 - u;

a = 2;
unom = @(x) x^2 - a;
x0 = 0;
fnc = @(t,x) dyn(x,unom(x));
[T,X] = ode45(fnc,[0 1], x0);

N = 100;
t = linspace(0,1,N);
Xt = interp1(T,X,t);
plot(t,a*t,'--',t,Xt,'-r');