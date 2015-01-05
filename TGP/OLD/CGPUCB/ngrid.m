function u = ngrid(ctr, a, b, size_r, size_theta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NONUNIFORM GRID
% create a nonuniform ellipsoidal grid using an exponential function
% on the radius to create a uniform mesh of delta(r)/r
%
% Inputs:
%
% ctr - center of the ellipse
% a,b - focal lengths of the ellipse
% size_r - number of discretizations of the radius
% size_theta - number of discretizations of the angle
%
% Output: control points sampled 
%
% Author : Okan Koc, ETHZ IDSC Lab, 2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ctr_x = ctr(1);
ctr_y = ctr(2);

r = linspace(0,1,size_r);
r = exp(r-1);
r = r - 1/exp(1);

u = zeros(2, size_r * size_theta);

for i = 1:size_r
    for j = 1:size_theta
        idx = (i-1) * size_theta + j;
        theta = j*2*pi / size_theta;
        u(1,idx) = ctr_x + a * r(i) * cos(theta);
        u(2,idx) = ctr_y + b * r(i) * sin(theta);
    end
end