% loop generation
r = 100;
h = 0.02;
N = 100;
t = linspace(h,h*N,N);
rr = r .* (t >= 0.5) .* (t <= 1.5);

z1 = r - r*cos(2*pi*t) + 2*rr .* (1 + cos(2*pi*t)); 
z2 = r*sin(2*pi*t);
%ddz1 = 4*pi*pi*(r-2*rr).*cos(2*pi*t);
%ddz2 = -4*pi*pi*r*sin(2*pi*t);

plot(z1,z2);