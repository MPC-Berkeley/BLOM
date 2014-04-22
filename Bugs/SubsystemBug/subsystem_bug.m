%% note: In this file, there seems to be a problem with the constraint of F_tilde at the first or last time step. 
% BLOM is not properly handling the constraint.
% run this file and it will show plots of where the constraints are
% violated


%% Problem Details

% load obstacle file here
%load('stopLightObstacles.mat')

% Initial Problem Formulation
N = 50;                % horizon length
x_0 = 0;
y_0 = 0;
phi_0 = 0;
v_0 = 0;
x_end = [0 50];             % final x position
y_end = [0 0];              % final y position
phi_end = [0 0];            % final phi position

% Constants for simulation
dt_model    =       0.2;              %   time step for model                           [s]

% load times when obstacles are on and off
% [stoplight_on] = generateStoplights(8,10,dt_model,300);
% times_on_actual = [ones(length(stoplight_on),2) repmat(stoplight_on,1,4)]; % generate array of times
% times_on_add_actual = -5*(times_on_actual==0);
% times_on_add_actual(:,1:2) = -5;

% Constants related to obstacles avoidance
dsafe       =       0.2;              %  safety distance obstacle avoidance             [m] 
K           =       3;                %  number of position states                      [1]
n0          =       size(obstacles,1);%  number of obstacles                            [1]



% Constants for car model
l           =       4.91;             %   vehicle length  (Hyundai grandeur specs)      [m]
l_r         =       3;                %   length rear end - cg                          [m]
l_f         =       l-l_r;            %   length front end - cg                         [m]
w           =       1.86;             %   width of vehicle (Hyundai grandeur specs)     [m]
t           =       0.15;             %   wheel thickness                               [m]
r           =       0.3;              %   wheel radius                                  [m] 
V_lb        =       -5;               %   lower bound for V                             [m/s] 
V_ub        =       13.4;             %   upper bound for V. 30 mph = 13.4 mph          [m/s]

beta_ub     =       (0.9)*(pi/2);     %   upper bound delta_f                           [rad]
beta_lb     =       -(0.9)*(pi/2);    %   lower bound delta_f                           [rad]
beta_r_ub   =       0.2*dt_model;     %   upper bound rate delta_f                      [rad/s]
beta_r_lb   =       -0.2*dt_model;    %   lower bound rate delta_f                      [rad/s]

F_tilde_ub  =       1500;             %   upper bound F_f                               [N]
F_tilde_lb  =       -1500;            %   lower bound F_f                               [N]
F_tilde_r_ub=       Inf;              %   upper bound rate F_f                          [N/s]
F_tilde_r_lb=       -Inf;             %   lower bound rate F_f                          [N/s]



beta_lb     =       0; % -pi/2*.6;         %   lower bound for beta                          [rad] 
beta_ub     =       0; %pi/2*.6;          %   upper bound for beta                          [rad] 
alpha_min   =       0; %-.2*dt_model;     %   maximum angular acceleration                  [rad/s] 
alpha_max   =       0; %.2*dt_model;      %   minimum angular acceleration                  [rad/s] 
a_min       =       -1*dt_model;      %   minimum acceleration                          [m/s^2]   
a_max       =       1*dt_model;       %   maximum acceleration                          [m/s^2]
I           =       3344;             %   vehicle inertia                               [kg m^2]
m_car       =       1920;             %   vehicle mass                                  [kg]



% constants for car model simulation (includes missmatch)
m_car_sim   =      m_car + 0.05*m_car;%  vehicle mass simulation model  (+200 kg)       [kg]
l_r_sim     =      l_r + 0.05*l_r;    %  vechicle rear end -cg simulation model (+ 5 cm)[m]                
l_f_sim     =      l_f + 0.05*l_f;    %  vechicle frontend-cg simulation model (+ 5 cm) [m]

% save results for plot_results
x_m                 =   [x_0,y_0,phi_0,0]';     % set first measurement as initial position
res_vec             =   x_m;
ui                  =   [0 0]';
input_vec           =   ui;
x_ol                =   zeros(1,N);
y_ol                =   zeros(1,N);
phi_ol              =   zeros(1,N);
V_ol                =   zeros(1,N);

% setup obstacles
make_robot_poly = @(x) orientedBoxToPolygon([x(1), x(2), l_f,l_r, w, x(3)]);


% for the first time step, do not consider obstacles, i.e. set the function
% g(x) to zero. 

x_jac       = zeros(1,n0+1);
y_jac       = zeros(1,n0+1);
phi_jac     = zeros(1,n0+1);
g_val_x0    = zeros(1,n0+1);
times_on    = zeros(1,n0+1);
times_on_add = -5*ones(1,n0+1);

% x_jac       = [0 0]; %zeros(1,n0+1);
% y_jac       = [0 0]; %zeros(1,n0+1);
% phi_jac     = [0 0]; %zeros(1,n0+1);
% g_val_x0    = [0 0]; %zeros(1,n0+1);

x_prev      = [0 0];
y_prev      = [0 0];
phi_prev    = [0 0];

% weights on cost function
weight_phi = [0 0];
weight_y = [0 0];
weight_x = [0 10; (N-3)*dt_model 0.5];
obstacle_penalty = [0 0];

% initial values 
F_tilde = [0 0; 1 0];
beta = [0 0; 0 0];
beta_init = 0; %initial condition for beta

%% Start Solving Optimization problem with BLOM

% load initial problem without obstacles to create an initial trajectory
current_model = 'Subsystem_ProjectModel';

BLOM_SetDataLogging(current_model)
ModelSpec = BLOM_ExtractModel(current_model,N,dt_model);
RunResults = BLOM_RunModel(ModelSpec);
SolverStruct = BLOM_ExportToSolver(ModelSpec,'IPOPT');
[OptGuess ExtVars InitialStates ] = BLOM_SplitResults(ModelSpec,RunResults);
SolverStructData =  BLOM_SetProblemData(SolverStruct,ModelSpec,OptGuess, ExtVars, InitialStates);
SolverResult  =  BLOM_RunSolver(SolverStructData,ModelSpec);

figure(1)
plot(SolverResult.x_dynamics)
figure(2)
plot(SolverResult.F_tilde)