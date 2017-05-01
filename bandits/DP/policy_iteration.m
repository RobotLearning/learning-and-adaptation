%% Policy Iteration for the infinite horizon case

function [V,pi] = policy_iteration(Pr,R,beta)

A = size(Pr,3);
n = size(Pr,1);
V = zeros(n,1);
Vcand = zeros(n,A);
pi = zeros(n,1);
pi_old = ones(n,1);
iter = 0;

while pi ~= pi_old
    iter = iter + 1;
    pi_old = pi;
    for i = 1:A
        Vcand(:,i) = sum(Pr(:,:,i) .* R(:,:,i),2) + beta * Pr(:,:,i) * V;
    end
    pi_old = pi;
    [~,pi] = max(Vcand,[],2);
    % iterate till convergence or solve linear equations for V
    V = (eye(n) + beta * Pr(:,:,pi)) \ sum(Pr(:,:,pi) .* R(:,:,pi),2);
    
end
disp(['Iterated ', num2str(iter), ' times!']);