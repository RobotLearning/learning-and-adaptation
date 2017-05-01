% Generate from an inverse gamma distribution
% parameterized using a,b

function xInvGamma = InvGam(a,b)

k = a;
theta = 1/b;

% First generate a gamma-distr. random var.
% Marsaglia's simple transformation-rejection method

d = k - 1/3;
c = 1/sqrt(9*d);
xUni = rand;
xNormal = randn;
v = (1 + c*xNormal)^3;
while v <= 0 || log(xUni) >= 0.5*xNormal^2 + d - d*v + d*log(v)
    xNormal = randn;
    xUni = rand;
    v = (1 + c*xNormal)^3;
end

xGamma = d*v;
% scale 
xGamma = theta*xGamma;
xInvGamma = 1/xGamma;


end