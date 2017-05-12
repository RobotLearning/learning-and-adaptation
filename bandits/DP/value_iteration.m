%% Value Iteration for the infinite horizon case

% V are the values and pi are the optimal policies
function [V,pi] = value_iteration(Pr,R,beta)

A = size(Pr,3);
n = size(Pr,1);
V = zeros(n,1);
delta_min = 1e-3;
delta = 1;
Q = zeros(n,A);
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

% Computes Q values
% If rewards are dependent on transition (s and s') or not (only s')
% uses Hadamard product and sums along second dimension
function Q = compute_q(P,R,V,beta) 

    n = size(P,1);
    A = size(P,3);
    Q = zeros(n,A);
    if size(R,3) > 1
        for i = 1:A
            Q(:,i) = sum(P(:,:,i) .* R(:,:,i),2) + beta * P(:,:,i) * V;
        end
    else
        for i = 1:A
            Q(:,i) = P(:,:,i) * (R(:,i) + beta * V);
        end
    end
end