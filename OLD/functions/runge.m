%% Investigating Runge's phenomenon

% Generate the Runge's function
runge_fnc = @(x) 1./(1 + 25*(x.^2));
n = 100;
t = linspace(-1,1,n+1);

% fit with Bernstein polynomials
n = 10;
x = linspace(-1,1,n+1);

% evaluate matrix A
% generate Bernstein bases
bernstein = @(x,v,n) (factorial(n)/(factorial(v)*factorial(n-v))) * ...
             (x.^v) .* ((1-x).^(n-v));         
A = zeros(n+1,n+1);
for i = 1:n+1
    A(:,i) = bernstein(x,i-1,n);
end
b = runge_fnc(x); 

% solve linear system
beta = A\b(:);

% evaluate interpolating polynomial
f = zeros(1,length(t));
for i = 1:n+1
    f = f + beta(i) * bernstein(t,i-1,n);
end

% spline
sp = spline(x,b,t);

plot(t,runge_fnc(t), '-r', t, f,'-b', t, sp, '-g');
legend('Original function', 'interpolating polynomial', 'cubic spline');

%% MATLAB Curve-fitting toolbox splines

n = 10;
x = linspace(0,1,n+1);
m = 100;
t = linspace(0,1,m+1);
y = sin(2*pi*x);

% derivatives set to 0 on endpoints
x2 = [x(1)*ones(1,4), x, x(end)*ones(1,4)];
y2 = [y(1), zeros(1,4), y(2:end), zeros(1,4)];

sp4 = spapi(4,x,y);
sp10 = spapi(10,x,y);
sp10_drv = spapi(10,x2,y2);

fnplt(sp4, '-r'); hold on;
pp10 = fn2fm(sp10,'pp');
f1 = ppval(pp10,t);
plot(t,f1,'-b'); hold on;
fnplt(sp10_drv, '--k');
plot(x,y,'o')
legend('cubic','9th order','9th order with drv', 'data')
hold off

%% Bezier curves 

n = 10;
x = linspace(0,1,n+1);
m = 100;
t = linspace(0,1,m+1);
y = sin(2*pi*x);
% derivative
dy = 2*pi*cos(2*pi*x);
% evaluate dy at endpoints
xx = [x; x];
x2 = reshape(xx,1,2*length(x));
dydy = [y; dy];
dy2 = reshape(dydy,1,2*length(x));
yy = spapi(3,x2,dy2);
pp = fn2fm(yy,'pp');
yy = ppval(pp,t);
% derivatives should be exact
plot(t,yy,'-r'); hold on

% bezier curve
% make up 10 bezier curves and then transform
% generate Bernstein bases
bstn = @(x,v,n) (factorial(n)/(factorial(v)*factorial(n-v))) * ...
             (x.^v) .* ((1-x).^(n-v)); 
         
bezier3 = @(x,p0,p1,dp0,dp1) p0*bstn(x,0,3)+(p0+dp0)*bstn(x,1,3)+...        
            (p1+dp1)*bstn(x,2,3)+p1*bstn(x,3,3);

yy_bez = [];
h = 1/n;
for i = 1:n
    vec = x(i):h:x(i+1);
    dp0 = dy(i)/(n*3);
    dp1 = dy(i+1)/(n*3);
    bez = bezier3(vec,x(i),x(i+1),dp0,dp1);
    yy_bez = [yy_bez bez(1:end-1)];
end
yy_bez(m+1) = y(end);
plot(t,yy_bez,'-g');
plot(x,y,'o')
legend('cubic spline with drv','bezier curve','data')
hold off