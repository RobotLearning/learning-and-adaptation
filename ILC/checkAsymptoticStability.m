% Function that checks for asymptotial stability of the nominal plant
% Based on the actual plant dynamics

function AS = checkAsymptoticStability(act,nom)

    F = nom;
    Z = act;
    
    E = Z - F;
    s = svd(F); % returns only singular values
    
    if norm(E) < min(s)
        AS = true;
        disp('Model is asymptotically stable!');
    else
        warning('Model may not be asymptotically stable!');
        AS = false;
    end


end