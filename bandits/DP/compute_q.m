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