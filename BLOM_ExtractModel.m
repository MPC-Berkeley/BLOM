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
% Output:
%   ModelSpec - BLOM optimization problem description
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

load_system(name); % TODO - replace it later to load_system, after got rid of gcs insid ExtractModel
[all_names, AAs ,  Cs , ineq ,cost ,in_vars,all_state_vars,ex_vars] = ExtractModel(horizon,dt,integ_method,name);

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

% merge everything into one array
ModelSpec.A = horzcat(cost.A', ineq.AAs{:})';
ModelSpec.ineq_start_A = size(cost.A,1)+1;
ModelSpec.ineq_end_A = size(ModelSpec.A,1);
ModelSpec.A = horzcat(ModelSpec.A', AAs{:})';
ModelSpec.eq_start_A = ModelSpec.ineq_end_A  +1;
ModelSpec.eq_end_A = size(ModelSpec.A,1);


ModelSpec.C = blkdiag(cost.C, ineq.Cs{:});
ModelSpec.ineq_start_C = size(cost.C,1)+1;
ModelSpec.ineq_end_C = size(ModelSpec.C,1);
ModelSpec.C = blkdiag(ModelSpec.C, Cs{:});
ModelSpec.eq_start_C = ModelSpec.ineq_end_C+1;
ModelSpec.eq_end_C = size(ModelSpec.C,1);

% do vectorization of all_names ahead of time
num_terms = cellfun(@length, strfind(all_names,';'))' + 1; % number of ';'
terms_so_far = [0; cumsum(num_terms)]; % is number of multiple names
all_fields = textscan([all_names{:}],'BL_%sOut%dt%d','Delimiter','.;');
vec_idx = zeros(terms_so_far(end),1); % preallocate vec_idx
vec_idx(terms_so_far(1:end-1)+1) = 1:length(all_names); % first of each
twoterms = find(num_terms == 2);
vec_idx(terms_so_far(twoterms)+2) = twoterms; % 2nd of each
multiterms = find(num_terms > 2); % multiples, should be fewer of these
for i = 1:length(multiterms)
    vec_idx(terms_so_far(multiterms(i))+2 : ...
        terms_so_far(multiterms(i)+1)) = multiterms(i);
end
ModelSpec.all_names_struct.terms_so_far = terms_so_far;
ModelSpec.all_names_struct.all_fields = all_fields;
ModelSpec.all_names_struct.vec_idx = vec_idx;
