function [OptGuess ExtVars InitialStates ] = BLOM_SplitResults(ModelSpec,RunResults)
%
% [OptGuess ExtVars InitialStates ] = BLOM_SplitResults(ModelSpec,RunResults)
%
% Splits full structure of results to 3 types of data: first guess
% variables, external parameters, initial states.
%
% Input:
%   ModelSpec   -   Model structure generatated by BLOM_ExtractModel.
%   RunResults  -   Structure containing variables data.
% 
% Output:
%   OptGuess    -   Structure with variables for an initial guess.
%   ExtVars     -   Structure with external variables.
%   InitialStates-  Structure with initial state variables.

ResultsVec = BLOM_ConvertStructToVector(ModelSpec.all_names,RunResults);

OptGuess =  BLOM_GetStructOfOptGuess(ModelSpec,ResultsVec);
ExtVars  = BLOM_GetStructOfExtVars(ModelSpec,ResultsVec);
InitialStates =  BLOM_GetStructOfInitialStates(ModelSpec,ResultsVec);