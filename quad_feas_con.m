function [c,ceq] = feas_con(z,PATH,PAR,CON)
%--------------------------------------------------------------------------
%This function evaluates all the nonlinear constraints at a finite number
%of time instances for the current choice of the optimization parameters.
%
%INPUTS
%   z:    decision variables: z := [T,l_points_init]
%           - T is the trajectory end time
%           - l_points are the INTERIOR lambda profile support points
%   PATH: struct containing the path spline and its derivatives wrt. LAMBDA
%           - PATH.p_spline: the path spline object
%           - PATH.dl_p_spline: first derivative wrt. LAMBDA
%           - PATH.ddl_p_spline: second derivative wrt. LAMBDA
%   PAR:  struct containing all system parameters
%   CON:  struct containing all constraint bounds
%
%OUTPUTS
%   c: vector that contains the nonlinear inequalities c(z)<= 0 evaluated
%      at given time instances 't'.
%   ceq: vector that contains the nonlinear equalities c(z) = 0 evaluated
%      at given time instances 't'.
%--------------------------------------------------------------------------

%% Read system and constraint parameters
a = 0.25;
b = PAR.Iy/(2*PAR.m*PAR.L);
g = PAR.g;

fmax          = CON.fmax;
fmin          = CON.fmin;
fi_max        = CON.fi_max;
fi_min        = CON.fi_min;
fi_dot_max    = CON.fi_dot_max;
phi_max      = CON.phi_max;
phi_dot_max  = CON.phi_dot_max;
phi_ddot_max = CON.phi_ddot_max;
NumT          = CON.NumT; %number of time points where constraints are evaluated

% p_spline   = PATH.p_spline;
dl_p_spline  = PATH.dl_p_spline;
ddl_p_spline = PATH.ddl_p_spline;
d3l_p_spline = PATH.d3l_p_spline;
d4l_p_spline = PATH.d4l_p_spline;
d5l_p_spline = PATH.d5l_p_spline;

%REMARK: Recall that the above derivatives are with respect to the path
%        parameter lambda and NOT with respect to time.
%        The path derivatives wrt TIME are computed below and are referred
%        to as dt_p_spline and ddt_p_spline.

%% Compute current lambda profile
%The lambda profile is computed depending on the current value of the
%optimization variables, that is, depending on 'T' (end time) and 
%'points_l', which is a vector containing the support points of the lambda
%spline.
%REMARK: The support points containted in 'p_middle' define the lambda
%        profile in the open interval (0,T), while the start point 
%        lambda(0) = 0 and the end point lambda(T) = 1 are fixed and, thus,
%        not subject to optimization.

%--------------------------------------------------------------------------
%Read current values of the optimization variables:
T = z(1); %this is the value of the current end time 

%--------------------------------------------------------------------------
%Build the lambda profile, given the values of the optimization variables

%The lambda profile support points, including start (0) and end (1) points.
l_points = [0,z(2:end),1];

%get equidistant time points on the interval [0,T] which provide the time 
%information corresponding to the support points in 'l_points'
sup_x = linspace(0,T,length(l_points));

%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%Compute the LAMBDA profile spline

%REMARK: zero first and second derivatives wrt. time of lambda are required
%        to guarantee zero first derivatives wrt. time of the y- and
%        z-paths.
%        This is done here by supplying 4 additional zeros after sup_y.
% sp = spapi(augknt(sup_x,6,1),...
%     [sup_x,sup_x(1),sup_x(1),sup_x(end),sup_x(end)],...
%     [l_points,0,0,0,0]);
sp = spapi(augknt(sup_x,10,1),...
    [sup_x,...
     repmat(sup_x(1),1,4),repmat(sup_x(end),1,4)],...
    [l_points,zeros(1,8)]);

l_spline = fn2fm(sp,'pp'); %convert to 'pp' form

%Differentiate
dt_l_spline  = fnder(l_spline,1); %first derivative wrt. time
ddt_l_spline = fnder(l_spline,2); %second derivative wrt. time
d3t_l_spline = fnder(l_spline,3); %third derivative wrt. time
d4t_l_spline = fnder(l_spline,4); %fourth derivative wrt. time
d5t_l_spline = fnder(l_spline,5); %fourth derivative wrt. time

%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%

%Discretize time:
%The constraints are evaluated at discrete points in time.
t_points = linspace(0,T,NumT);

% figure(21)
% plot(t_points,ppval(l_spline,t_points),'-b',xl,points_lambda,'or')
% title('Initial Lambda Profile: from inside optimization')
% xlabel('real time [s]')
% ylabel('lambda')

%% Evaluate splines

%--------------------------------------------------------------------------
%First we evaluate the lambda spline at the TIME instances given by the
%vector 't_points'.
lambda     = ppval(l_spline,t_points);
dt_lambda  = ppval(dt_l_spline,t_points);
ddt_lambda = ppval(ddt_l_spline,t_points);
d3t_lambda = ppval(d3t_l_spline,t_points);
d4t_lambda = ppval(d4t_l_spline,t_points);
d5t_lambda = ppval(d5t_l_spline,t_points);

%--------------------------------------------------------------------------
%Now we evaluate the path splines at the LAMBDA VALUES!!!
%REMARK: Path splines are always evaluated at lambda values, never at time
%        values.
%REMARK: The following '*_path' are matrices of dimension (2, length(t_points)).
%        The first row corresponds to the y values, the second row to the z
%        values.

dl_path  = ppval(dl_p_spline, lambda);
ddl_path = ppval(ddl_p_spline, lambda);
d3l_path = ppval(d3l_p_spline, lambda);
d4l_path = ppval(d4l_p_spline, lambda); 
d5l_path = ppval(d5l_p_spline, lambda);
%the prefix ddl means 'second derivative with respect to lambda'

%--------------------------------------------------------------------------
%Compute the TIME DERIVATIVES of the path
%(according to the chain rule)
%REMARK: dt_lambda and ddt_lambda have size(1,length(t_points)). Thus, we
%        replicate these vectors once vertically.

%second derivative wrt. TIME of the path spline:
ddt_path = repmat(dt_lambda.^2,2,1).*ddl_path + ...
             dl_path.*repmat(ddt_lambda,2,1);
         
%third derivative wrt. TIME of the path spline
d3t_path = 3.*repmat(dt_lambda,2,1).*ddl_path.*repmat(ddt_lambda,2,1) + ...
               repmat(dt_lambda.^3,2,1).*d3l_path + ...
               dl_path.*repmat(d3t_lambda,2,1);

%fourth derivative wrt. TIME of the path spline
d4t_path = 6.*repmat(dt_lambda.^2,2,1).*repmat(ddt_lambda,2,1).*...
               d3l_path + ddl_path.*(3.*repmat(ddt_lambda.^2,2,1) +...
               4.*repmat(dt_lambda,2,1).*repmat(d3t_lambda,2,1)) + ...
               repmat(dt_lambda.^4,2,1).*d4l_path + ...
               dl_path.*repmat(d4t_lambda,2,1);
           
%fifth derivative wrt. TIME of the path spline
%fifth derivative wrt. TIME of the path spline
d5t_path = 10.*ddl_path.*repmat(ddt_lambda,2,1).*repmat(d3t_lambda,2,1) + ...
               10.*repmat(dt_lambda.^2,2,1).*d3l_path.*repmat(d3t_lambda,2,1) + ...
               10.*repmat(dt_lambda.^3,2,1).*repmat(ddt_lambda,2,1).*d4l_path + ...
               5.*repmat(dt_lambda,2,1).*(3.*repmat(ddt_lambda.^2,2,1).* ...
               d3l_path + ddl_path.*repmat(d4t_lambda,2,1)) + ...
               repmat(dt_lambda.^5,2,1).*d5l_path + ...
               dl_path.*repmat(d5t_lambda,2,1);
           
%rename
ddt_y = ddt_path(1,:); %second derivative wrt. TIME of y-path evaluated at 
                       %the TIME points 't_points'.
ddt_z = ddt_path(2,:); %second derivative wrt. TIME of z-path evaluated at
                       %the TIME points 't_points'.
            

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INEQUALITY CONSTRAINT FUNCTIONS    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RECALL: the nonlinear constraint functaion have to be of the form
%          c(x)  <= 0

%% Constraints on common thrust
%--------------------------------------------------------------------------
%Compute the common thrust
fc = sqrt(ddt_y.^2 + (ddt_z + g).^2);

%Compute first derivative wrt. TIME of common thrust
dt_fc = (ddt_y.*d3t_path(1,:)+(PAR.g + ddt_z).*...
             d3t_path(2,:)) ./ fc;
         
%Compute second derivative wrt. TIME of common thrust
ddt_fc = (-(2.*ddt_y.*d3t_path(1,:) + 2.*(PAR.g + ...
             ddt_z).*d3t_path(2,:)).^2 + ...
             4.*(ddt_y.^2 + (PAR.g+ddt_z).^2).*...
             (d3t_path(1,:).^2 + d3t_path(2,:).^2 + ddt_y.*...
             d4t_path(1,:) + (PAR.g + ddt_z).*d4t_path(2,:)))./...
             (4.*(ddt_y.^2 + (PAR.g+ddt_z).^2).^(3/2) );
         
%Compute third derivative wrt. TIME of common thrust
d3t_fc = (3.*(2.*ddt_path(1,:).*d3t_path(1,:) + ...
             2.*(PAR.g + ddt_path(2,:)).*d3t_path(2,:)).^3 - ...
             12.*(ddt_path(1,:).^2 + (PAR.g + ddt_path(2,:)).^2).* ...
             (2.*ddt_path(1,:).*d3t_path(1,:) + 2.*(PAR.g + ddt_path(2,:)).* ...
             d3t_path(2,:)).*(d3t_path(1,:).^2 + d3t_path(2,:).^2 + ...
             ddt_path(1,:).*d4t_path(1,:) + (PAR.g + ddt_path(2,:)).* ...
             d4t_path(2,:)) + ...
             4.*(ddt_path(1,:).^2 + (PAR.g + ddt_path(2,:)).^2).^2.* ...
             (6.*d3t_path(1,:).*d4t_path(1,:) + ...
             6.*d3t_path(2,:).*d4t_path(2,:) + ...
             2.*ddt_path(1,:).*d5t_path(1,:) + 2.*(PAR.g + ddt_path(2,:)).*...
             d5t_path(2,:)))./ ...
             (8.*(ddt_path(1,:).^2 + (PAR.g + ddt_path(2,:)).^2).^(5/2));
         
%--------------------------------------------------------------------------
%Constraints on common thrusts 
%(fc <= f_max) <--> (fc - f_max <= 0)
c_fmax = fc - fmax;

%(fc >= f_min) <--> (-fc <= -f_min) <--> (-fc + fmin <= 0)
c_fmin = -fc + fmin;


%% Constraints on anglular velocity and angular acceleration
%--------------------------------------------------------------------------
%Compute the roll angle
beta = asin(ddt_y./fc);

%Compute first derivative wrt. TIME of beta
dt_beta = (-dt_fc.*ddt_y + fc.*d3t_path(1,:))./...
              (fc.^2.*sqrt(1-(ddt_y.^2)./(fc.^2)));
          
%Compute second derivative wrt. TIME of beta
ddt_beta = (-dt_fc.^2.*ddt_y.^3 + fc.*ddt_fc.*...
                ddt_y.^3 - fc.^3.*(ddt_fc.*ddt_y+...
                2.*dt_fc.*d3t_path(1,:)) + fc.^4.*d4t_path(1,:)+...
                fc.^2.*ddt_y.*(2.*dt_fc.^2 + d3t_path(1,:).^2 - ...
                ddt_y.*d4t_path(1,:)) )./...
                (fc.^3.*(fc.^2-ddt_y.^2).*...
                sqrt(1-ddt_y.^2./(fc.^2)) );
            
%Compute third derivative wrt. TIME of beta
d3t_beta = (-2.*dt_fc.^3.*ddt_path(1,:).^5 + ...
               3.*fc.*dt_fc.*ddt_fc.*ddt_path(1,:).^5 + ...
               fc.^2.*(5.*dt_fc.^3.*ddt_path(1,:).^3 - ...
               ddt_path(1,:).^5.*d3t_fc) - ...
               fc.^6.*(ddt_path(1,:).*d3t_fc + ...
               3.*ddt_fc.*d3t_path(1,:) + 3.*dt_fc.*d4t_path(1,:)) + ...
               fc.^4.*ddt_path(1,:).*(-6.*dt_fc.^3 + ddt_path(1,:).* ...
               (2.*ddt_path(1,:).*d3t_fc + 3.*ddt_fc.*d3t_path(1,:)) + ...
               dt_fc.*(-9.*d3t_path(1,:).^2 + 3.*ddt_path(1,:).*d4t_path(1,:))) + ...
               fc.^7.*d5t_path(1,:) + ...
               fc.^5.*(6.*dt_fc.*ddt_fc.*ddt_path(1,:) + ...
               6.*dt_fc.^2.*d3t_path(1,:) + d3t_path(1,:).^3 + ...
               3.*ddt_path(1,:).*d3t_path(1,:).*d4t_path(1,:) - ...
               2.*ddt_path(1,:).^2.*d5t_path(1,:)) + ...
               fc.^3.*ddt_path(1,:).^2.* ...
               (-9.*dt_fc.*ddt_fc.*ddt_path(1,:) + ...
               3.*dt_fc.^2.*d3t_path(1,:) + 2.*d3t_path(1,:).^3 - ...
               3.*ddt_path(1,:).*d3t_path(1,:).*d4t_path(1,:) + ...
               ddt_path(1,:).^2.*d5t_path(1,:))) ./ ...
               (fc.^4.*(fc.^2 - ddt_path(1,:).^2).^2 .* ...
               sqrt(1 - ddt_path(1,:).^2./(fc.^2)));
 
%--------------------------------------------------------------------------
%Constraints on angle
%(beta <= beta_max) -> (beta - beta_max <= 0)
c_bmax = beta - phi_max;

%(beta >= -beta_max) -> (-beta <= beta_max) -> (-beta - beta_max <= 0)
c_bmin = -beta - phi_max;
           
%--------------------------------------------------------------------------
%Constraints on angular velocity (corresponds to second input u(2))
%(beta_dot <= beta_dot_max) --> (beta_dot - beta_dot_max <= 0)
c_dbmax = dt_beta - phi_dot_max;
 
%(beta_dot >= -beta_dot_max) --> (-beta_dot <= beta_dot_max) --> (-beta_dot - beta_dot_max <= 0)
c_dbmin = -dt_beta - phi_dot_max;

%--------------------------------------------------------------------------
%Constraint on angular acceleration 
%(corresponds to the second input's rate of change u_dot(2))

%(ddt_beta <= beta_ddot_max) --> (ddt_beta - beta_ddot_max <= 0)
c_ddbmax = ddt_beta - phi_ddot_max;

%(ddt_beta >= -beta_ddot_max) --> (-ddt_beta <= beta_ddot_max) --> (-ddt_beta - beta_ddot_max <= 0)
c_ddbmin = -ddt_beta - phi_ddot_max;

 
%% Stricter input constraints
%--------------------------------------------------------------------------
%Compute single thrusts
f1 = a.*fc + b.*ddt_beta;
f3 = a.*fc - b.*ddt_beta;

%Compute first derivative
dt_f1 = a.*dt_fc + b.*d3t_beta;
dt_f3 = a.*dt_fc - b.*d3t_beta;
%--------------------------------------------------------------------------
%Constraints on single thrusts
%(f1 <= fi_max) --> (f1 - fi_max <= 0)
c_f1max = f1 - fi_max;

%(f1 >= fi_min) --> (-f1 <= -fi_min) --> (-f1 + fi_min <= 0)
c_f1min = -f1 + fi_min;

%(f3 <= fi_max) --> (f3 - fi_max <= 0)
c_f3max = f3 - fi_max;

%(f3 >= fi_min) --> (-f3 <= -fi_min) --> (-f3 + fi_min <= 0)
c_f3min = -f3 + fi_min;

%--------------------------------------------------------------------------
%Constraints on single thrusts rate of change
%(dt_f1 <= fi_dot_max) --> (dt_f1 - fi_dot_max <= 0)
c_dt_f1max = dt_f1 - fi_dot_max;

%(dt_f1 >= -fi_dot_max) --> (-dt_f1 <= fi_dot_max) -->
%(-dt_f1 - fi_dot_max <= 0)
c_dt_f1min = -dt_f1 - fi_dot_max;

%(dt_f3 <= fi_dot_max) --> (dt_f3 - fi_dot_max <= 0)
c_dt_f3max = dt_f3 - fi_dot_max;

%(dt_f3 >= -fi_dot_max) --> (-dt_f3 <= fi_dot_max) -->
%(-dt_f3 - fi_dot_max <= 0)
c_dt_f3min = -dt_f3 - fi_dot_max;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   RETURN VECTORS                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c = [c_f1max, c_f1min, ...
     c_f3max, c_f3min, ...
     c_dt_f1max, c_dt_f1min, ...
     c_dt_f3max, c_dt_f3min, ...
     c_bmax, c_bmin, ...
     c_dbmax, c_dbmin, ...
     c_ddbmax, c_ddbmin, ...
     -dt_lambda,...
%      -ddt_lambda
     ];

ceq = []; %no equality constraints

end %EOF

