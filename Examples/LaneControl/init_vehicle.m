%% Control Parameters
dt = 0.1;      % Sampling time (control)
% dt_ol = 0.001;
% DT = 0.05;     % Sampling time (numerical integration)
DT_plant = 0.005;
% DT_anal = 0.05;
nz = 6;         % No. of states
nu = 3;         % No. of inputs
Hu = 15;        % Control horizon
Hi = 1;         % No. of time steps over which input is kept constant    
Hp = Hu*Hi;     % Prediction horizon

%% Vehicle Parameters
a = 1.432;          % CoG front axle distance [m]
b = 1.472;          % CoG rear axle distance [m]
c = 0.8125;         % Half trackwidth of vehicle [m]
r = 0.335;          % Wheel radius [m]
m = 2050;           % Vehicle mass [kg]
Iz = 3344;          % Moment of inertia about vertical axis [kg-m^2]
Iw_f = 1.2;         % Moment of inertia of wheel [kg-m^2]
Iw_r = 1.2;         % Moment of inertia of wheel [kg-m^2]
mu = 1;             % Road friction coefficient
cof = mu;
g = 9.81;           % Acceleration due to gravity [m/s^2]
C_f = 6.5e4;        % Front tire cornering stiffness [N]  
C_r = 6.5e4;        % Rear tire cornering stiffness [N]
Fz_f = b*m*g/(a+b)/2;  % Front tire static normal force [N]
Fz_r = a*m*g/(a+b)/2;  % Rear tire static normal force [N]
Fz = [Fz_f, Fz_f, Fz_r, Fz_r];
delta_f_min = -0.5; % Minimum front steering angle [rad]
delta_f_max = -delta_f_min;  % Maximum front steering angle [rad]
ddelta_max = 60*pi/180;
Tb_min = -0.5*Iw_f*mu*g/r;       % Minimum braking torque at each wheel [N-m]
Tb_max = -Tb_min;        % Maximum braking torque at each wheel [N-m]
Fl_min = -mu*m*g/4; % Minimum longitudinal force [N]
Fl_max = -Fl_min;   % Maximum longitudinal force [N]
beta_min = -0.5;    % Minimum braking ratio
beta_max = 0.5;     % Maximum braking ratio
alpha_lim = 10*pi/180;  % Maximum slip angle

%for mu=0.5
% C_lin      = -4.4448e+004;
% % C_sat      = 185.3451;
% max_alpha  = 0.05;

%for mu=1
% C_lin      = -5.3458e+004;
% max_alpha  = 0.08;
% offset_sat = -max_alpha*(C_lin-C_sat);
