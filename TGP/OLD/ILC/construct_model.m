function MATS = construct_model(t,x,PAR,u_trj,CON)

dim = size(x,1);
N = size(x,2);
Nu = N - 1;
h = t(2) - t(1);
dim_u = size(u_trj,1);

%%%%%%%%%%%% CONSTRUCT A and B matrices %%%%%%%%%%%%%%%%%%%%%%%
A = zeros(dim,dim,N);
B = zeros(dim,dim_u,N);
for i = 1:N
    [~,A(:,:,i), B(:,:,i)] = quadrocopterNominalDynamics(t(i),x(:,i),PAR,u_trj(:,i),true);
    % get discrete approximation from jacobian
    % crude approximation
    %A(:,:,i) = eye(dim,dim) + h * A(:,:,i);
    %B(:,i) = h * B(:,i);
    % exact matrix calculation 
    Mat = [A(:,:,i), B(:,:,i); zeros(dim_u, dim + dim_u)];
    MD = expm(h * Mat);
    A(:,:,i) = MD(1:dim,1:dim);
    B(:,:,i) = MD(1:dim,dim+1:end);
end

%%%%%%%%% CONSTRUCT lifted domain matrices F, G, and H %%%%%%%%
% G is identity matrix, H is zero matrix
F = zeros(N*dim, Nu*dim_u); % u is two-dimensional
G = eye(N*dim); % since y(t) = x(t) 
H = zeros(N*dim,Nu*dim_u); 
for l = 1:N
    for m = 1:Nu
        vec_x = (l-1)*dim + 1: l*dim;
        vec_u = (m-1)*dim_u + 1:m*dim_u;
        if m <= l
            MATS = eye(dim);
            for i = m+1:l
                MATS = MATS * A(:,:,i);
            end
            F(vec_x,vec_u) = MATS * B(:,:,m); % on diagonals only B(:,m)
        end
    end
end

% lifted domain constraints

%%%%%%%%%%%%%%%%%%%% CONSTRUCT umin and umax %%%%%%%%%%%%%
% extract constraints

umin(1,:) = CON.fmin - u_trj(1,1:Nu);
umin(2,:) = -CON.phi_dot_max - u_trj(2,1:Nu);
umax(1,:) = CON.fmax - u_trj(1,1:Nu);
umax(2,:) = CON.phi_dot_max - u_trj(2,1:Nu);

% arrange them in a format suitable for optimization
umin = umin(:);
umax = umax(:);

%%%%%%%%%%% CONSTRUCT L AND q %%%%%%%%%%%%%%
% construct D
D = (diag(ones(1,dim_u*(Nu-1)),dim_u) - eye(dim_u*Nu))/h;
D = D(1:end-dim_u,:); % D is (Nu-1)*nu x Nu*nu dimensional
% construct L1 and L2
L1 = zeros(Nu-1, Nu*dim_u);
L2 = zeros(Nu-1, Nu*dim_u);
a = 1/4;
b = PAR.Iy/(2*PAR.m*PAR.L); 
b = b/h; %b_bar
vec1 = [a -b 0 b];
vec2 = [a b 0 -b];
for i = 1:Nu-1
    L1(i,:) = [zeros(1,(i-1)*dim_u), vec1, zeros(1,(Nu-i-1)*dim_u)];
    L2(i,:) = [zeros(1,(i-1)*dim_u), vec2, zeros(1,(Nu-i-1)*dim_u)];
end
u_dot_max = [4*CON.fi_dot_max; CON.phi_ddot_max];
U_dot_max = repmat(u_dot_max,Nu-1,1);
u_star = u_trj(:); u_star = u_star(1:end-dim_u);

L = [D; -D; L1; -L1; L2; -L2];
q = [U_dot_max - D*u_star; 
     U_dot_max + D*u_star;
     CON.fmax - L1*u_star;
     -CON.fmin + L1*u_star;
     CON.fmax - L2*u_star;
     -CON.fmin + L2*u_star;];     

% return structure
MATS.F = F; 
MATS.G = G; 
MATS.H = H;
MATS.umin = umin; 
MATS.umax = umax; 
MATS.L = L; 
MATS.q = q;