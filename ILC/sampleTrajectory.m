% Sample a trajectory in phase space
% drawn from a multivariate normal distribution
% with squared exponential kernel based covariance matrix
%
% we condition on x0, the initial state

function traj = sampleTrajectory(n,N,x0)

mu = 0;
x = (1:N)/N;
l = 0.2;
s = 1.0;
kern = @(x1,x2) s * exp(-(norm(x1(:)-x2(:),2)^2)/(2*(l^2)));
Sigma = zeros(N);
for i = 1:N
    for j = 1:N
        Sigma(i,j) = kern(x(i),x(j));
    end
end

% condition on 0 at 0
Sigma = condition(Sigma,kern);

% spectral decomposition instead 
[V,D] = eig(Sigma);
D(D < 0) = 0.0;
% cholesky decomposition does not work for smooth kernels
%noise = 0;
%traj = mu + chol(Sigma + noise*eye(N)) * randn(N,1);
traj = mu + V*sqrt(D)*randn(N,n);

end

% Conditioning the Gaussian Process on y at x0 
% TODO: include more general case
function Sigma = condition(Sigma,kern)

    N = size(Sigma,1);
    x = (1:N)/N;
    kern_vec = zeros(1,N);
    for i = 1:N
        kern_vec(i) = kern(0,x(i));
    end
    Sigma = Sigma - kern(0,0)*(kern_vec'*kern_vec);
end