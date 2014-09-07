function [t,u_trj,x] = traj_gen(trj,PAR,h)

N = length(trj);
t = linspace(0,1,N);
A = PAR.A;
B = PAR.B;
dim_x = length(A);
dim = size(B,2);

Mat = [A, B; zeros(dim, dim_x + dim)];
MD = expm(h * Mat);
A_d = MD(1:dim_x,1:dim_x);
B_d = MD(1:dim_x,dim_x+1:end);

% nominal trajectory to be followed
x = zeros(dim_x,N);
x(2,1) = trj(2,1);
u_trj = zeros(dim,N);

% minor TODO: this can also be done without a for loop in matlab
for i = 1:N-1
    % deviation along desired trajectory
    delta = trj(2,i+1) - A_d(2,:)*trj(:,i);
    % solve for u
    u_trj(i) = delta / B_d(2,1);
    % simulate trj
    %x(:,i+1) = A_d * x(:,i) + B_d * u_trj(i);
end

u_trj(:,N) = u_trj(:,N-1);
% simulate nominal trajectory with RK4
handle = @springDamperNominalDynamics;
x = sim_RK4(t,x(:,1),u_trj,0,PAR,handle);

end