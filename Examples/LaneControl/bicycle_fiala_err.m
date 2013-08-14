function [statesdot, varargout] = bicycle_fiala_err(z,u,curv,m,mu,Iz,a,b,c,Fz_f,Fz_r,C_f,C_r)
%% Continuous time nonlinear bicycle model in error coordinates
% x: states:[xdot, ydot, Psidot, e_psi, e_y, s]'
% u: input:[delta_f, beta_fl, beta_fr, beta_rl, beta_rr]'
% curv = d(psi)/ds

% States
xdot = z(1);
ydot = z(2);
psid = z(3);
e_psi = z(4);
e_y = z(5);
s = z(6);

% Inputs
delta_f = u(1);
delta_r = 0;
beta_fl = u(2);
beta_fr = u(3);
beta_rl = u(4);
beta_rr = u(5);
beta = u(2:5);
eta = sqrt(1-beta.^2);

% Normal forces
Fz = [Fz_f; Fz_f; Fz_r; Fz_r];

% Cornering stiffnesses
C_alpha = [C_f; C_f; C_r; C_r];

% sin/cos
cos_f = cos(delta_f);
sin_f = sin(delta_f);
cos_e_psi = cos(e_psi);
sin_e_psi = sin(e_psi);

% Wheel speeds
vy_f = ydot+a*psid;
vy_r = ydot-b*psid;
vx_l = xdot-c*psid;
vx_r = xdot+c*psid;

vc_fl = vy_f*cos_f-vx_l*sin_f;
vc_fr = vy_f*cos_f-vx_r*sin_f;
vc_rl = vy_r*1-vx_l*0;
vc_rr = vy_r*1-vx_r*0;

vl_fl = vy_f*sin_f+vx_l*cos_f;
vl_fr = vy_f*sin_f+vx_r*cos_f;
vl_rl = vy_r*0+vx_l*1;
vl_rr = vy_r*0+vx_r*1;

% Slip angles
alpha_fl = atan(vc_fl/vl_fl);
alpha_fr = atan(vc_fr/vl_fr);
alpha_rl = atan(vc_rl/vl_rl);
alpha_rr = atan(vc_rr/vl_rr);
alpha = [alpha_fl; alpha_fr; alpha_rl; alpha_rr];

% Return slip angles if needed
if nargout > 1
    varargout{1} = [alpha_fl; alpha_fr; alpha_rl; alpha_rr];
end

% Longitudinal forces
Fl = mu*beta.*Fz;
Fl_fl = Fl(1);
Fl_fr = Fl(2);
Fl_rl = Fl(3);
Fl_rr = Fl(4);

% Limiting slip angles
alpha_sl = atan(3*mu*eta.*Fz./C_alpha);
if nargout > 2
    varargout{2} = alpha_sl;
end

% Cornering forces
Fc = zeros(length(beta),1);
for ii=1:length(Fc)
    if abs(alpha(ii)) < alpha_sl(ii)
        Fc(ii) = -C_alpha(ii)*tan(alpha(ii)) + C_alpha(ii)^2/(3*eta(ii)*mu*Fz(ii))*abs(tan(alpha(ii)))*tan(alpha(ii)) - C_alpha(ii)^3/(27*eta(ii)^2*mu^2*Fz(ii)^2)*tan(alpha(ii))^3;
    else
        Fc(ii) = -eta(ii)*mu*Fz(ii)*sign(alpha(ii));
    end
end
Fc_fl = Fc(1);
Fc_fr = Fc(2);
Fc_rl = Fc(3);
Fc_rr = Fc(4);
if nargout > 3
    varargout{3} = Fc;
end

% Body axes forces
Fx_fl = Fl_fl*cos_f - Fc_fl*sin_f;
Fy_fl = Fl_fl*sin_f + Fc_fl*cos_f;
Fx_fr = Fl_fr*cos_f - Fc_fr*sin_f;
Fy_fr = Fl_fr*sin_f + Fc_fr*cos_f; 
Fx_rl = Fl_rl*1 - Fc_rl*0;
Fy_rl = Fl_rl*0 + Fc_rl*1;
Fx_rr = Fl_rr*1 - Fc_rr*0;
Fy_rr = Fl_rr*0 + Fc_rr*1; 

% State derivatives
dsdt = (xdot*cos_e_psi - ydot*sin_e_psi)*(1/(1-e_y*curv));
statesdot(1,1) = psid*ydot + (1/m)*(Fx_fl + Fx_fr + Fx_rl + Fx_rr);
statesdot(2,1) = -psid*xdot + (1/m)*(Fy_fl + Fy_fr + Fy_rl + Fy_rr);
statesdot(3,1) = (1/Iz)*(a*(Fy_fl+Fy_fr) - b*(Fy_rl+Fy_rr) + c*(-Fx_fl + Fx_fr - Fx_rl + Fx_rr));
statesdot(4,1) = psid - curv*dsdt;
statesdot(5,1) = xdot*sin_e_psi + ydot*cos_e_psi;
statesdot(6,1) = dsdt;