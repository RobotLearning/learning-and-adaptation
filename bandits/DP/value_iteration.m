%% Value Iteration for the infinite horizon case

% V are the values and pi are the optimal policies
function [V,pi] = value_iteration(Pr,R,beta)

A = size(Pr,3);
n = size(Pr,1);
V = zeros(n,1);
delta_min = 1e-3;
delta = 1;
Vcand = zeros(n,A);
pi = zeros(n,1);

iter = 0;
while delta > delta_min
    
    for i = 1:A
        Vcand(:,i) = sum(Pr(:,:,i) .* R(:,:,i),2) + beta * Pr(:,:,i) * V;
    end
    [Vnew,pi] = max(Vcand,[],2);
    delta = max(abs(Vnew - V));
    V = Vnew;
    iter = iter + 1;
end
disp(['Iterated ', num2str(iter), ' times!']);