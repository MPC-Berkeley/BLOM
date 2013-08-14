function zdot = bicycle_model_linear(t,z)

% States
x_dot = z(1);
y_dot = z(2);
psi_dot = z(3);
e_psi = z(4);
e_y = z(5);
s = z(6);

% Inputs
delta_f = -0.01*e_y-0.1*e_psi;
beta_l = 0;
beta_r = 0;

% Derivatives (model obtained from BLOM)
zdot = zeros(6,1);
% zdot(1) = y_dot*psi_dot+4.905*beta_l+4.905*beta_r;
% zdot(2) = (-x_dot*x_dot*psi_dot+126.8293*y_dot-2.5366*psi_dot-63.4146*x_dot*delta_f)/(x_dot); 
% zdot(3) = (-1.555*y_dot+163.9543*psi_dot-55.6699*x_dot*delta_f-2.4431*x_dot*beta_l+2.4431*x_dot*beta_r)/(x_dot);
% zdot(4) = psi_dot;
% zdot(5) = x_dot*sin(e_psi)+y_dot*cos(e_psi);
% zdot(6) = x_dot*cos(e_psi)-y_dot*sin(e_psi);

% Derivatives (model from scratch)
init_vehicle
alpha_f = (y_dot+a*psi_dot)/x_dot-delta_f;
alpha_r = (y_dot-b*psi_dot)/x_dot;
zdot(1) = y_dot*psi_dot + (cof/m)*(beta_l*(Fz(1)+Fz(3)) + beta_r*(Fz(2)+Fz(4)));
zdot(2) = -x_dot*psi_dot + (1/m)*(-2*C_f*alpha_f - 2*C_r*alpha_r);
zdot(3) = (1/Iz)*(-2*a*C_f*alpha_f + 2*b*C_r*alpha_r) + (c*cof/Iz)*(-beta_l*(Fz(1)+Fz(3))+beta_r*(Fz(2)+Fz(4)));
zdot(4) = psi_dot;
zdot(5) = x_dot*sin(e_psi) + y_dot*cos(e_psi);
zdot(6) = x_dot*cos(e_psi) - y_dot*sin(e_psi);