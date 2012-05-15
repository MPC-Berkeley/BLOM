function ModelSpec = BLOM_ExtractModel(name,horizon,dt,integ_method,options)
% ModelSpec = BLOM_ExportModel(name,horizon,dt,disc,options)
%
% Converts simulink model to BLOM optimization problem.
%
% Input: 
%   name    -   name of the model to convert
%   horizon -   integer number of steps in prediction horizon length, 1 if not supplied
%   dt      -   time step size [s], 1 if not supplied  
%   disc    -   discretization method {'none','Euler','Trapez','RK4'}
%   options -   options created by BLOM_optset function.
%
%Examples:
%   ModelSpec = BLOM_ExportModel; 
%       Converts current active model.
%
%   ModelSpec = BLOM_ExportModel('Sys',10,0.5);  
%       Converts current model 'Sys', with a a time step of 0.5 sec and 10 time steps.
%
%   ModelSpec = BLOM_ExportModel('Sys',10,0.5,);  
%       Converts current model 'Sys', with a a time step of 0.5 sec and 10 time steps, 
%       assuming discrete model.
%
%   ModelSpec = BLOM_ExportModel('Sys',10,0.5,'Trapez');  
%       Converts current model 'Sys', with a a time step of 0.5 sec ,10
%       time steps and trapezoidal discretization method.


if (nargin < 1)
    name = gcs;
end
if(nargin < 2)
    horizon = 1;
end
if(nargin < 3)
    dt = 1;
end
if(nargin < 4)
    integ_method = 'none';
end
if(nargin < 5)
    options = [];
end

open_system(name); % TODO - replace it later to load_system, after got rid of gcs insid ExtractModel
[all_names, AAs ,  Cs , ineq ,cost ,in_vars,all_state_vars,ex_vars] = ExtractModel(horizon,dt,integ_method);

ModelSpec.name = name;
ModelSpec.integ_method = integ_method;
ModelSpec.dt = dt;
ModelSpec.horizon = horizon;
ModelSpec.options = options;
ModelSpec.all_names = all_names;
ModelSpec.AAs = AAs;
ModelSpec.Cs = Cs;
ModelSpec.ineq = ineq;
ModelSpec.cost = cost;
ModelSpec.in_vars = in_vars;
ModelSpec.all_state_vars = all_state_vars;
ModelSpec.ex_vars = ex_vars;

