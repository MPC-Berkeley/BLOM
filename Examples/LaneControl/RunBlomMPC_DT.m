%% MPC for obstacle avoidance using BLOM
clear
% close all
clc

%% Vehicle parameters
init_vehicle

%% Road parameters
road_len = 0:0.1:500;
road_wid = 5;
curv = 0;

%% Object data
obj(1).s = 40;
obj(1).e_y = -0.5;
obj(1).off_s = 10;
obj(1).off_e_y = 1.5;

%% Simulation parameters
simLength = 50;        % Total simulation time steps
Zinit = [15 0 0 0 0 0]';     % [xdot ydot Psidot e_psi e_y s];
Zcurr = Zinit;
Zlb = [0,-10,-10,-pi/2,-road_wid,0]';
Zub = [25,10,10,pi/2,road_wid,road_len(end)]';
Zref = [15,0,0,0,0,0];
Ulb = [delta_f_min, beta_min*ones(1,2)]';
Uub = [delta_f_max, beta_max*ones(1,2)]';
Q = [1 1 1 0 20 0];
R = [5 5 5];

%% BLOM parameters
% Vehicle model state derivatives
zdot_f_A = [0 1 1 zeros(1,6);...    % ydot*psidot
            zeros(1,7) 1 0;...      % beta_l
            zeros(1,8) 1;...        % beta_r
            2 0 1 zeros(1,6);...    % xdot^2*psidot
            0 1 zeros(1,7);...      % ydot
            0 0 1 zeros(1,6);...    % psidot
            1 zeros(1,5) 1 0 0;...  % xdot*delta
            1 zeros(1,6) 1 0;...    % xdot*beta_l
            1 zeros(1,7) 1;...      % xdot*beta_r
            0 0 1 0 1 zeros(1,4);...    % psidot*ey
            1 0 0 BLOM_FunctionCode('cos') zeros(1,5);...   % xdot*cos(epsi)
            0 1 0 BLOM_FunctionCode('sin') zeros(1,5);...   % ydot*sin(epsi)
            1 0 0 BLOM_FunctionCode('sin') zeros(1,5);...   % xdot*sin(epsi)
            0 1 0 BLOM_FunctionCode('cos') zeros(1,5);...   % ydot*cos(epsi)
            ];
zdot_f_C = [1, mu*(Fz(1)+Fz(3))/m, mu*(Fz(2)+Fz(4))/m, zeros(1,11);...
            zeros(1,3), -1, (-2*C_f-2*C_r)/m, (-2*C_f*a+2*C_r*b)/m, 2*C_f/m, zeros(1,7);...
            zeros(1,4), (-2*a*C_f+2*b*C_r)/Iz, (-2*a^2*C_f-2*b^2*C_r)/Iz, 2*a*C_f/Iz, -c*mu*(Fz(1)+Fz(3))/Iz, c*mu*(Fz(2)+Fz(4))/Iz, zeros(1,5);...
            zeros(1,5), 1, zeros(1,3), -curv, -curv, curv, zeros(1,2);...
            zeros(1,12), 1, 1;...
            zeros(1,10), 1, -1, 0, 0;...
            ];
zdot_g_A = [zeros(1,9);...  % constant
            1 zeros(1,8);...  % xdot
            zeros(1,4) 1 zeros(1,4);... % ey
            ];
zdot_g_C = [1 0 0; 0 1 0; 0 1 0; 1 0 -curv; 1 0 0; 1 0 -curv];

% Input constraints
Uconstr_f_A = [eye(nu); zeros(1,nu)];
Uconstr_f_C = [eye(nu), -Uub; -eye(nu), Ulb];
Uconstr_g_A = zeros(1,nu);
Uconstr_g_C = 1;

% State constraints
Zconstr_f_A = [eye(nz);...          % states 
               zeros(1,nz);...      % constant
               zeros(1,5), 2;...    % s^2
               zeros(1,4), 2, 0;... % e_y^2
               ];
Zconstr_f_C = [eye(nz), -Zub, zeros(nz,2);...   % state upper bounds
               -eye(nz), Zlb, zeros(nz,2);...   % state lower bounds
               zeros(1,4), 2*obj(1).e_y/obj(1).off_e_y^2, 2*obj(1).s/obj(1).off_s^2,...
                        1 - obj(1).s^2/obj(1).off_s^2 - obj(1).e_y^2/obj(1).off_e_y^2, -1/obj(1).off_s^2, -1/obj(1).off_e_y^2;...  % obstacle avoidance
               ];
Zconstr_g_A = zeros(1,nz);
Zconstr_g_C = 1;

% Front slip angle constraint
Sconstr_f_A = [1, zeros(1,8);...            % xdot
               0, 1, zeros(1,7);...         % ydot
               0, 0, 1, zeros(1,6);...      % psidot
               1, zeros(1,5), 1, 0, 0;...   % delta*xdot
               ];
Sconstr_f_C = [-alpha_lim, 1, a, -1;...
               -alpha_lim, -1, -a, 1];
Sconstr_g_A = [1, zeros(1,8)];  % xdot
Sconstr_g_C = [1; 1];

% Cost
cost_f_A = [2*eye(nz+nu); eye(nz), zeros(nz,nu)];
cost_f_C = [Q, R, -2*Q.*Zref];
cost_g_A = zeros(1,nz+nu);
cost_g_C = 1;

%% BLOM
ModelSpec = BLOM_ExtractModel('BLOM_MPC_DT',Hp/dt);
RunResults = BLOM_RunModel(ModelSpec);
SolverStruct = BLOM_ExportToSolver(ModelSpec,'IPOPT');
[OptGuess, ExtVars, InitialStates] = BLOM_SplitResults(ModelSpec,RunResults);

%% MPC
Z_all = zeros(nz,simLength);
dt_ipopt = -1*ones(1,simLength);
for n=1:simLength
    
    %% Set initial guess, initial states
    Z_all(:,n) = Zcurr;
    if n==1
%         OptGuess.System_u = zeros(Hp,nu);
%         OptGuess.System_z = ones(Hp,1)*Zcurr';
%         OptGuess.System_zdot = zeros(Hp,nz);
    else
%         OptGuess.System_u = SolverResult.System_u;
%         OptGuess.System_z = SolverResult.System_z;
%         OptGuess.System_zdot = SolverResult.System_zdot;
%         OptGuess.System_slip_angle_constraint = SolverResult.System_slip_angle_constraint;
        OptGuess = SolverResult;
    end
%     OptGuess.System_u = zeros(Hp,nu);
%     OptGuess.System_z = ones(Hp,1)*Zcurr';
%     OptGuess.System_zdot = zeros(Hp,nz);
    InitialStates.System_z = Zcurr';
    
    %% Run IPOPT
    SolverStructData =  BLOM_SetProblemData(SolverStruct,ModelSpec,OptGuess,ExtVars,InitialStates);
    options.print_level = 0;
    options.print_user_options = 'no';
    tic
    [SolverResult, ResultsVec, ResultInfo] = BLOM_RunSolver(SolverStructData,ModelSpec,options);
    dt_ipopt(n) = toc;

    %% Get optimal solution
    Uopt = SolverResult.System_u(1,:)';
    Ufiala = [Uopt(1), Uopt(2), Uopt(3), Uopt(2), Uopt(3)]';
    
    %% Simulate plant
    Zcurr = bicycle_fiala_err_discrete(Zcurr,Ufiala,curv,dt,DT_plant,m,mu,Iz,a,b,c,Fz_f,Fz_r,C_f,C_r);    
    
    %% Plot results
    figure(10)
    subplot(211)
    plot(SolverResult.System_z(:,1)); %hold on
    ylabel('$$\dot{x} [m/s]$$','interpreter','latex')
    xlabel('Time step')
    if n==simLength
        hold off
    end
    
    subplot(212)
    plot(SolverResult.System_z(:,6),SolverResult.System_z(:,5),'b','linewidth',2); hold on;
    plot(SolverResult.System_z(:,6),road_wid*ones(Hp,1),'r');
    plot(SolverResult.System_z(:,6),-road_wid*ones(Hp,1),'r');
    plot(SolverResult.System_z(:,6),zeros(Hp,1),'--b');
    ell_x = (obj(1).s-obj(1).off_s):0.01:(obj(1).s+obj(1).off_s);
    ell_y1 = obj(1).e_y + obj(1).off_e_y*sqrt(1 - (ell_x-obj(1).s).^2/obj(1).off_s^2);
    ell_y2 = obj(1).e_y - obj(1).off_e_y*sqrt(1 - (ell_x-obj(1).s).^2/obj(1).off_s^2);
    plot(ell_x,ell_y1,'r','linewidth',2)
    plot(ell_x,ell_y2,'r','linewidth',2)
    hold off
    if n==simLength
        hold off
    end
    ylabel('$$Y [m]$$','interpreter','latex');
    xlabel('$$X [m]$$','interpreter','latex')
    
    figure(11)
    subplot(211);
    plot(SolverResult.System_u(:,1)); %hold on
    ylabel('$$\delta [rad]$$','interpreter','latex');
    if n==simLength
        hold off
    end
    
    subplot(212)
    plot(SolverResult.System_u(:,2)); hold on
    plot(SolverResult.System_u(:,3));
    hold off
    ylabel('$$\beta$$','interpreter','latex')
    xlabel('Time step')
    if n==simLength
        hold off
    end    
end

%% Figures
figure(1)
subplot(211)
plot(0:dt:dt*(simLength-1),Z_all(1,:),'b','linewidth',2)
ylabel('$$\dot{x} [m/s]$$','interpreter','latex')
xlabel('$$Time [s]$$','interpreter','latex')

subplot(212)
plot(Z_all(6,:),Z_all(5,:),'b','linewidth',2); hold on;
plot(Z_all(6,:),road_wid*ones(simLength,1),'r');
plot(Z_all(6,:),-road_wid*ones(simLength,1),'r');
plot(Z_all(6,:),zeros(simLength,1),'--b');
ell_x = (obj(1).s-obj(1).off_s):0.01:(obj(1).s+obj(1).off_s);
ell_y1 = obj(1).e_y + obj(1).off_e_y*sqrt(1 - (ell_x-obj(1).s).^2/obj(1).off_s^2);
ell_y2 = obj(1).e_y - obj(1).off_e_y*sqrt(1 - (ell_x-obj(1).s).^2/obj(1).off_s^2);
plot(ell_x,ell_y1,'r','linewidth',2)
plot(ell_x,ell_y2,'r','linewidth',2)
ylabel('$$Y [m]$$','interpreter','latex');
xlabel('$$X [m]$$','interpreter','latex')

