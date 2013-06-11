function [RunResults ResultsVec]= BLOM_RunModel(ModelSpec,options)
%
%  [RunResults ResultsVec]= BLOM_RunModel(ModelSpec,options)
%
%   Executes the simulink model and returns the recorded values.
%   NOTE: Changes base workspace variables. 
%
% Input:
%   ModelSpec -   Model structure generatated by BLOM_ExtractModel.
%   options   -   options created by BLOM_optset function.
%
% Output:
%   RunResults -    Structure with fields according to ModelSpec, holding
%                   the simulation results. 
%   ResultsVec -    Vector with the same results    

sim(ModelSpec.name,0:ModelSpec.dt:ModelSpec.dt*ModelSpec.horizon);
logsout.unpack('all');

%num_terms = cellfun(@length, strfind(ModelSpec.all_names,';')) + 1;
%terms_so_far = [0, cumsum(num_terms)];
%all_fields = textscan([ModelSpec.all_names{:}],'BL_%sOut%dt%d','Delimiter','.;');
terms_so_far = ModelSpec.all_names_struct.terms_so_far;
all_fields = ModelSpec.all_names_struct.all_fields;
base_name = all_fields{1}(terms_so_far(1:end-1) + 1); % only first name
output_number = all_fields{2}(terms_so_far(1:end-1) + 1); % only first name
time_index = all_fields{3}(terms_so_far(1:end-1) + 1); % only first name
port_number = all_fields{4}(terms_so_far(1:end-1) + 1); % only first name
vec_index = all_fields{5}(terms_so_far(1:end-1) + 1); % only first name


varName = cell(size(base_name));
for ii = 1:length(varName)
    varName{ii} = [base_name{ii} '_' num2str(port_number(ii))];
end

% convert from Simulink ToWS to BLOM vector one variable name at a time
% goal is to copy each matrix of (time, port) data in a vectorized way
[sorted, index] = sort(varName);
% First sort by names, then find last occurrence of each name
[names, I] = unique(sorted);
I = [0; I];

ResultsVec = zeros(length(ModelSpec.all_names),1);
% very ugly : everything is done in the base workspace
for i=1:length(names)
    % look for the variable in the base workspace
    if (isempty(who(names{i})))
        warning(['Var ' names{i} ' not found in workspace']);
        continue;
    end
    data_i = eval([names{i} '.Data']);
    data_i(1,:) = [];
    
    % convert times and port #'s for this signal name into 1d indices
    time_indices_i = time_index(index(I(i)+1:I(i+1)));
    vec_indices_i = vec_index(index(I(i)+1:I(i+1)));
    inds_i = sub2ind(size(data_i), time_indices_i, vec_indices_i);
    ResultsVec(index(I(i)+1:I(i+1))) = data_i(inds_i);
end
RunResults = BLOM_ConvertVectorToStruct(ModelSpec.all_names_struct,ResultsVec);

%{
ResultsVec1 = zeros(length(ModelSpec.all_names),1);
for i=1:length(ModelSpec.all_names)
    vname = ModelSpec.all_names{i};
    vname = vname(1:min([length(vname), strfind(vname,';')-1]));
    idx = strfind(vname,'.');
    if length(idx)~=2
        continue;
    end
    name = vname(1:idx(1)-1);
    % look for the variable in the base workspace
    if (isempty(evalin('base',['who(''' name ''')'])))
        warning(['Var ' name ' not found in base workspace']);
        continue;
    end
    % Time index 
    time = vname(idx(2)+2:end);
    % Variable index for vector variables
    port = vname(idx(1)+4:idx(2)-1);
    
    %wsname = sprintf('%s.signals.values(%s,%s)', name, time, port);
    wsname = [name '.signals.values(' time ',' port ')'];
    % Take the variable from the base workspace
    ResultsVec1(i) = evalin('base', wsname);
end
if ~isequal(ResultsVec, ResultsVec1)
    disp('mismatch')
end
%}
