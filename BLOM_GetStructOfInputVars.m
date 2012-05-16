function  InpVars = BLOM_GetStructOfInputVars(ModelSpec,ResultsVec)
%
% EmptyExtVars = BLOM_GetStructOfExtVars(ModelSpec)
%
% Creates template structure with the externally specified parameters.
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


if (nnz(ModelSpec.in_vars) == 0)
    InpVars = struct;
else
    vec = ones(nnz(ModelSpec.in_vars ~= 0 ),1); 
     if (nargin < 2)
        ResultsVec = nan*vec;
    else
       ResultsVec = ResultsVec(ModelSpec.in_vars ~= 0);
    end
    % TODO : the select vector should point to the original external
    % variable, traced back to ExtractModel
    InpVars = BLOM_ConvertVectorToStruct(ModelSpec.all_names(ModelSpec.in_vars ~= 0 ),ResultsVec,vec);
end
