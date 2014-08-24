function bounds = update_bnds(~,STR,u_last)

%--------------------------------------------------------------------------
% Updates the bounds for the optimization problem in CGP-UCB. 
% Uses derivative constraints to update the bounds from last iteration
%
%--------------------------------------------------------------------------

% extract components from structure
u = u_last; % used for 1-derivative constraints
h = STR.h;
CON = STR.CON;
fmin = CON.fmin;
fmax = CON.fmax;
f_dot_max = CON.f_dot_max;
phi_dot_max = CON.phi_dot_max;
phi_ddot_max = CON.phi_ddot_max;

% build bounds
f_last = u(1);
phi_dot_last = u(2);
f_next_min = f_last - h * f_dot_max;
f_next_max = f_last + h * f_dot_max;
phi_next_min = phi_dot_last - h * phi_ddot_max;
phi_next_max = phi_dot_last + h * phi_ddot_max;

% make sure strict bound constraints are also enforced
low_bnd = max([fmin,       -phi_dot_max; 
               f_next_min, phi_next_min]);
upp_bnd = min([fmax,       phi_dot_max; 
               f_next_max, phi_next_max]);

bounds = [low_bnd(:), upp_bnd(:)];
