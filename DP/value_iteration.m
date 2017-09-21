%% Value Iteration for the infinite horizon case

% V are the values and pi are the optimal policies
function [V,pi] = value_iteration(Pr,R,beta)

n = size(Pr,1);
V = zeros(n,1);
delta_min = 1e-3;
delta = 1;
pi = zeros(n,1);

iter = 0;
while delta > delta_min
    
    Q = compute_q(Pr,R,V,beta);        
    [Vnew,pi] = max(Q,[],2);
    delta = max(abs(Vnew - V));
    V = Vnew;
    iter = iter + 1;
end
disp(['Iterated ', num2str(iter), ' times!']);

end