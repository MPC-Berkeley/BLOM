function  EmptyOptGuess = BLOM_GetStructOfOptGuess(ModelSpec,ResultsVec)
%
% EmptyOptGuess = BLOM_GetStructOfInitialStates(ModelSpec)
%
% Creates template structure with the initial guess for the solver.
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

vec = ones(size(ModelSpec.all_names));
% TODO : the select vector should point to the original external
% variable, traced back to ExtractModel
if (nargin < 2)
    ResultsVec = nan*vec;
end

EmptyOptGuess = BLOM_ConvertVectorToStruct(ModelSpec.all_names_struct,ResultsVec,vec);
