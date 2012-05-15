function  EmptyInitialStates = BLOM_GetStructOfInitialStates(ModelSpec,ResultsVec)
%
% EmptyInitialStates = BLOM_GetStructOfInitialStates(ModelSpec)
%
% Creates template structure with the initial states of the system.
% This structure should be later filled by the user.
%
% Input:
%   ModelSpec -   Model structure generatated by BLOM_ExtractModel.
%   ResultsVec-   An optional vector with data, same size as ModelSpec.all_names
%
% Output:
%   EmptyExtVars -  Template structure with the externally specified
%                   parameters. Filled with NaNs if ResultsVec is not
%                   present


if (nnz(ModelSpec.all_state_vars) == 0)
    EmptyInitialStates = struct;
else
    vec = ones(nnz(ModelSpec.all_state_vars ~= 0 ),1); 
    if (nargin < 2)
        ResultsVec = nan*vec;
    else
       ResultsVec = ResultsVec(ModelSpec.all_state_vars ~= 0);
    end

    % TODO : the select vector should point to the original external
    % variable, traced back to ExtractModel
    EmptyInitialStates = BLOM_ConvertVectorToStruct(ModelSpec.all_names(ModelSpec.all_state_vars ~= 0 ),ResultsVec,vec);
end
