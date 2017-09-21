%% Policy Iteration for the infinite horizon case

function [V,pi] = policy_iteration(Pr,R,beta)

n = size(Pr,1);
pi = zeros(n,1);
pi_old = ones(n,1);
iter = 0;
V = zeros(n,1);

while pi ~= pi_old
    
    iter = iter + 1;
    pi_old = pi;    
    
    % COMPUTE POLICY
    Q = compute_q(Pr,R,V,beta);
    [~,pi] = max(Q,[],2);
    
    % UPDATE VALUES
    % rearrange to solve linear equation
    V = compute_v(Pr,R,pi,beta);
end
disp(['Iterated ', num2str(iter), ' times!']);

end

% Compute Values corresponding to updated policy
function V = compute_v(P,R,pi,beta)
    
    n = size(P,1);
    PP = zeros(n);
    % If rewards are dependent on transition (s and s') or not (only s')
    if size(R,3) > 1
        RR = zeros(n);
        for s = 1:n
            PP(s,:) = P(s,:,pi(s));
            RR(s,:) = R(s,:,pi(s));
        end

        % iterate till convergence or solve linear equations for V
        V = (eye(n) - beta * PP) \ sum(PP .* RR,2);
    else
        RR = zeros(n,1);
        for s = 1:n
            PP(s,:) = P(s,:,pi(s));
            RR(s) = PP(s,:) * R(:,pi(s));
        end
        V = (eye(n) - beta * PP) \ RR;
    end
end