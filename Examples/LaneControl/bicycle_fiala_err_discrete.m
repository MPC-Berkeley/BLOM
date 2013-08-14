function z1 = bicycle_fiala_err_discrete(z,u,curv,dt,DT,m,mu,Iz,a,b,c,Fz_f,Fz_r,C_f,C_r)
%% Discretize NL differential equation
t=0;
while t < dt
    current_dt = min(DT, dt-t);
    z = z + current_dt*bicycle_fiala_err(z,u,curv,m,mu,Iz,a,b,c,Fz_f,Fz_r,C_f,C_r);
    t = t+current_dt;
end
z1 = z;