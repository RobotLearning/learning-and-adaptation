function [t,u_trj,x] = quad_traj_gen(shape,PAR,CON,h,varargin)

% 2D Trajectory Generation with Splines: Main File
%
% 29.04.2011
% Fabian Mueller
% Modification date : 23.11.2012 and 16.12.2013
% by: OKAN KOC - simplified input and output, commented out plots
%
% INPUTS:
% 	shape :  string specifying the type of trajectory generated
% 	PAR  :   parameters specifying system constants
% 	CON  :   constraints on system dynamics
%   h    :   discretization constant
%   vargin : required y-z coordinates can be specified as well

%% Generate y and z path geometries

%--------------------------------------------------------------------------
% Procedure:
%   1) Choose an arbitrary number of points in the yz-plane. Write the
%      coordinates into the vectors 'y_coord' and 'z_coord'.
%   2) In the 'index' vector, specify the points where you want to define
%      the first and second spline derivative.
%   3) For the points specified by 'index', provide the derivative values
%      in the vectors 'dy' and 'dz'
%--------------------------------------------------------------------------

%-----------------
% Some example trajectories
%-----------------
switch shape
    case 'S shape'
        y_coord = [0, 0.3, -0.3, 0.1];
        z_coord = 1.5*[0, 0.3, 0.6, 0.9];        
        index = [1,4];
        dy = [0,0,1,0];
        dz = [0,1,0,1];
    case 'Horizontal' 
        % Point to point horizontal
        y_coord = [0, 1.5];
        z_coord = [0, 0];        
        index = [1,2];
        dy = [0,0,0,0]; %[dy(0),ddy(0),dy(1),ddy(1)]
        dz = [0,0,0,0]; %[dz(0),ddz(0),dz(1),ddz(1)]
    case 'Diagonal' 
        % Point to point diagonal
        y_coord = [0 2];
        z_coord = [0, 1.5];        
        index = [1,2];
        dy = [0,0,0,0]; %[dy(0),ddy(0),dy(1),ddy(1)]
        dz = [0,0,0,0]; %[dz(0),ddz(0),dz(1),ddz(1)]
    case 'Arc'
        y_coord = [0.2 1 2];
        z_coord = [0, 1.5, 0.5];        
        index = [1,3];
        dy = [0,0,0,0]; %[dy(0),ddy(0),dy(1),ddy(1)]
        dz = [0,0,0,0]; %[dz(0),ddz(0),dz(1),ddz(1)]
    case 'Wave'
        y_coord = [0 1 2 3];
        z_coord = [0, 1.5, 0 1.5];        
        index = [1,4];
        dy = [0,0,0,0]; %[dy(0),ddy(0),dy(1),ddy(1)]
        dz = [0,0,0,0]; %[dz(0),ddz(0),dz(1),ddz(1)]
    case 'Circle'
        index = [1,5];
        r=1;
        z_coord = r*[0,-1,0,1,0];
        y_coord = r*[1,0,-1,0,1];        
        dz = r*[-2*pi,0,-2*pi,0];
        dy = r*[0,-(2*pi)^2,0,-(2*pi)^2];
    case 'Question mark'
        y_coord = 1.8*[0, 0, -0.3, 0, 0.3];
        z_coord = 2*[0, 0.35, 0.65, 0.9, 0.65];        
        index = [1,5];
        dy = [0,0,0,0];
        dz = [1,1,-1,0];
    case 'A shape'
        y_coord = 1.8*[0, 0.2, 0.6, 0.57, 0.3,  -0.1];
        z_coord = 1.8*[0, 0.8, 0.1, 0.22, 0.4,  0.5];        
        index = [1,6];
        dy = [0,0,0,1];
        dz = [1,1,0,1];
    otherwise
        error('Shape not recognized. Please input a valid shape...');
end

% if y and z coordinates are specified extra, then extract them
if length(varargin) == 2
    y_coord = varargin{1};
    z_coord = varargin{2};
end

%% Compute the path splines

%Compute the splines that interpolate the given path points
y   = [y_coord;z_coord];
drv = [dy;dz];
lambda = linspace(0,1,size(y,2)); %equidistant support points
pd = [lambda(index(1)),lambda(index(2))];

%Compute 5th order spline
sp = spapi(augknt(lambda,6,1),[lambda,augknt(pd,2)],[y,drv]);

%Convert spline to 'pp' form
p_spline = fn2fm(sp,'pp');

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Evaluate the path splines
yy = ppval(p_spline, linspace(0,1,101)); %evaluate spline in y-z plane
dyy = ppval(fnder(p_spline,1), linspace(0,1,101));
ddyy = ppval(fnder(p_spline,2), linspace(0,1,101));

%Compute analytical derivatives of y,z wrt lambda
PATH.p_spline = p_spline;
PATH.dl_p_spline = fnder(p_spline,1);  %first analytical derivative 
PATH.ddl_p_spline = fnder(p_spline,2); %second analytical derivative
PATH.d3l_p_spline = fnder(p_spline,3); %third analytical derivative
PATH.d4l_p_spline = fnder(p_spline,4); %fourth analytical derivative
PATH.d5l_p_spline = fnder(p_spline,5); %fourth analytical derivative

%% Generate Initial Lambda Profile
%--------------------------------------------------------------------------
% Procedure
%   1) Provide an initial guess for the trajectory's end time: T_init
%   2) Define the number of interior lambda profile support points: N
%      (5 has proven to be a good choice)
%--------------------------------------------------------------------------

T_init = 2.5; %initial guess for trajectory's end time
N = 5;        %number of INTERIOR lambda profile support points
              %REMARK: N also equals the number of optimization variables

% Get the initial lambda profile and the support points
% 3rd order profile to provide initial points
% In a first step we find a minimum-order spline that goes through the
% points (0,0) and (T_init,1). 2 points -> minimum order is 3.
xx = [0,T_init];
yy = [0,1];

%Find a spline that interpolates the points [xx,yy]
lambda_3 = spline(xx,[0,yy,0]);

%Evaluate this spline at N+2 equidistant support points
sup_x = linspace(0,T_init,N+2); %equidistant support point time
sup_y = ppval(lambda_3,sup_x);  %the support point values

xl_init = sup_x; % initial support points (used for plot)
l_points_init = sup_y(2:end-1); %only consider interior points

%Get the time derivatives at the interior support points
drv1 = ppval(fnder(lambda_3,1),sup_x(2:end-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute a 9th order spline that interpolates the support points given by
%[sup_x,sup_y] and that has zero {first,second,third,fourth}-order
%derivative at the start and end point.
sp = spapi(augknt(sup_x,10,2),[sup_x,...
     repmat(sup_x(1),1,4),repmat(sup_x(end),1,4),sup_x(2:end-1)],...
    [sup_y,zeros(1,8),drv1]);

%Conversion from 'sp' to 'pp'
% l_spline is a spline in 'pp' form describing the lambda profile
l_spline = fn2fm(sp,'pp');
% l_points_init are the corresponding support points

%% Solve Optimization Problem
%Initial guess for decision variables
z0 = [T_init, l_points_init];

options = optimset('Display','off','MaxIter',100,'MaxFunEvals',1000, ...
                    'Algorithm', 'active-set');

%Minimize 'objfun' with 'z0' as the initial guess for the decision
%variables subject to the constraints defined in 'mycon'.
%We simply want to minimize the end time, which corresponds to the first
%element of the optimization variables vector z
%
%INPUTS
%   z: vector of optimization variables
%OUTPUTS
%   f: objective value, corresponds to the trajectory's end time T
%--------------------------------------------------------------------------
objfun = @(z) z(1);
[z,fval,exitflag] = fmincon(objfun,z0,[],[],[],[],[],[],...
                            @(z)quad_feas_con(z,PATH,PAR,CON),options);
%%%%%%%%%%%%%%%%%%%%%%                        
% Read optimization results:
T_opt = z(1);                  %optimized trajectory end time
l_points_opt = [0,z(2:end),1]; %the optimal lambda support points
supt_opt = linspace(0,T_opt,length(l_points_opt)); %corresponding time points

%Build the 9th-order lambda spline from the optimal support points
sp = spapi(augknt(supt_opt,10,1),[supt_opt,...
     repmat(supt_opt(1),1,4),repmat(supt_opt(end),1,4)],...
    [l_points_opt,zeros(1,8)]);

l_spline_opt = fn2fm(sp,'pp'); %convert to 'pp' form

%%%%%%%%%%%%%%%%%%%%%%
%Evaluate spline at discrete time instances
t_eval_opt = 0:h:T_opt;

lambda_opt     = ppval(l_spline_opt,t_eval_opt); %optimal lambda
dt_lambda_opt  = ppval(fnder(l_spline_opt,1),t_eval_opt); %first time derivative
ddt_lambda_opt = ppval(fnder(l_spline_opt,2),t_eval_opt); %second time derivative
d3t_lambda_opt = ppval(fnder(l_spline_opt,3),t_eval_opt); %third time derivative
d4t_lambda_opt = ppval(fnder(l_spline_opt,4),t_eval_opt); %fourth time derivative
d5t_lambda_opt = ppval(fnder(l_spline_opt,5),t_eval_opt); %fifth time derivative

%% Plot the initial and final lambda profiles

%Evaluate the lambda profile
t_eval_init = linspace(0,T_init,100);
lambda_init = ppval(l_spline,t_eval_init);

% figure
% h1 = plot(t_eval_init,lambda_init,'LineWidth',1,'LineStyle','--','Color','k'); hold all
% plot(xl_init,[0,l_points_init,1],'ok');
% h2 = plot(t_eval_opt,lambda_opt,'LineWidth',1,'LineStyle','-','Color','k');
% h3= plot(supt_opt,l_points_opt,'ok');
% hold off; grid on;
% title('Initial and final Lambda Profiles')
% legend([h1,h2,h3],'Initial lambda profile','Final lambda profile','Support points',4);
% xlabel('time [s]')
% ylabel('\lambda')

%% Compute states and thrust

%Evaluate the path along the initial and the optimized lambda profile
path_init = ppval(p_spline,lambda_init);

%derivative wrt TIME of lambda profile
path_opt     = ppval(p_spline,lambda_opt);
dl_path_opt  = ppval(fnder(p_spline,1),lambda_opt);
ddl_path_opt = ppval(fnder(p_spline,2),lambda_opt);
d3l_path_opt = ppval(fnder(p_spline,3),lambda_opt);
d4l_path_opt = ppval(fnder(p_spline,4),lambda_opt);
d5l_path_opt = ppval(fnder(p_spline,5),lambda_opt);

%first derivative wrt. TIME of the path spline
dt_path_opt = dl_path_opt.*repmat(dt_lambda_opt,2,1);

%second derivative wrt. TIME of the path spline
ddt_path_opt = repmat(dt_lambda_opt.^2,2,1).*ddl_path_opt + ...
               dl_path_opt.*repmat(ddt_lambda_opt,2,1);
           
%third derivative wrt. TIME of the path spline
d3t_path_opt = 3.*repmat(dt_lambda_opt,2,1).*ddl_path_opt.*repmat(ddt_lambda_opt,2,1) + ...
               repmat(dt_lambda_opt.^3,2,1).*d3l_path_opt + ...
               dl_path_opt.*repmat(d3t_lambda_opt,2,1);
           
%--------------------------------------------------------------------------
%Compute common thrust
fc_opt = sqrt(ddt_path_opt(1,:).^2+(ddt_path_opt(2,:)+PAR.g).^2);
     
%Compute first derivative wrt. TIME of common thrust
dt_fc_opt = (ddt_path_opt(1,:).*d3t_path_opt(1,:)+(PAR.g + ddt_path_opt(2,:)).*...
             d3t_path_opt(2,:)) ./ fc_opt;
         
%--------------------------------------------------------------------------
%Compute the roll angle
phi_opt = -asin(ddt_path_opt(1,:)./fc_opt);

%Compute first derivative wrt. TIME of beta
dt_phi_opt = (dt_fc_opt.*ddt_path_opt(1,:) - fc_opt.*d3t_path_opt(1,:) )./...
              (fc_opt.^2.*sqrt( 1-(ddt_path_opt(1,:).^2)./(fc_opt.^2) ) );

%% Plot y and z trajectories

% figure
% subplot(3,2,1); 
% plot(t_eval_opt,path_opt(1,:),'-r'); hold on; grid on;
% xlabel('time [s]')
% ylabel('y value [m]')
% title('State trajectory y in real time')
% legend('nominal y trajectory');
% 
% subplot(3,2,2); 
% plot(t_eval_opt,path_opt(2,:),'-r'); hold on; grid on;
% xlabel('time [s]')
% ylabel('z value [m]')
% title('State trajectory z in real time')
% legend('nominal z trajectory');

% Plot path in yz-plane
% subplot(3,2,[3 4]); 
% plot(path_opt(1,:),path_opt(2,:),'-or'); hold on; grid on;
% plot(y_coord,z_coord,'og');
% xlabel('y')
% ylabel('z')
% title('Path in y-z plane');
% legend('nominal y-z path');
% axis equal

% Plot angle trajectory
% subplot(3,2,3)
% h1 = plot(t_eval_opt,phi_opt); grid on; hold all;
% h2 = plot([t_eval_opt(1),t_eval_opt(end)],[CON.phi_max,CON.phi_max],'--r');
% plot([t_eval_opt(1),t_eval_opt(end)],[-CON.phi_max,-CON.phi_max],'--r');
% xlabel('time [s]')
% ylabel('Angle [rad]')
% title('Roll angle and angle rate trajectories')
% legend([h1,h2],'Angle [rad]','Constraints')

%% Plot control input : common thrust trajectory and roll rate

% subplot(2,1,1)
% fc = plot(t_eval_opt,fc_opt, '-r'); hold on; grid on;
% title('Collective thrust')
% xlabel('time [s]')
% ylabel('thrust [m/s^2]')
% 
% subplot(2,1,2)
% rollrate = plot(t_eval_opt,dt_phi_opt, '-r'); hold on; grid on;
% title('Roll rate')
% xlabel('time [s]')
% ylabel('Angular velocity [rad/s]')

%% Prepare States, Control Inputs and Time Optimized

t = t_eval_opt;
u_trj = [fc_opt; dt_phi_opt];
x = [path_opt(1,:); dt_path_opt(1,:); path_opt(2,:); dt_path_opt(2,:); phi_opt];
