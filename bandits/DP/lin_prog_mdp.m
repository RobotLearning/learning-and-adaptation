%% Linear program to solve infinite-horizon MDP

function [V,pi] = lin_prog_mdp(Pr,R,beta)

a = size(Pr,3);
n = size(Pr,1);
% prob distribution for V
mu = 1/n * ones(n,1);

P = reshape(permute(Pr,[1,3,2]),n*a,n);
I = repmat(eye(n),a,1);
A = beta * P - I;

b = zeros(n,a);
for i = 1:a
    b(:,i) = -Pr(:,:,i) * R(:,i);
end
b = b(:);

opt = optimoptions('linprog','Algorithm','Dual-simplex');
V = linprog(mu,A,b,[],[],[],[],opt);

Q = compute_q(Pr,R,V,beta);
[~,pi] = max(Q,[],2);