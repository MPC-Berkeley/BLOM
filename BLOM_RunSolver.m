function [ SolverResult  ResultsVec ]=  BLOM_RunSolver(SolverStruct,ModelSpec,options)
%
%   [ SolverResult  ResultsVec ]=
%   BLOM_RunSolver(SolverStruct,ModelSpec,options)
%
%   Runs the solver and returns the solution .
%
% Input:
%   SolverStruct - Solver description struct, created with BLOM_ExportToSolver.
%   ModelSpec -    Model structure generatated by BLOM_ExtractModel.
%   options   -    options created by BLOM_optset function.
%
% Output:
%   SolverResult -  Structure with fields according to ModelSpec, holding
%                   the solver results. 
%   ResultsVec -    Vector with the same results    


switch (SolverStruct.solver)
    case 'fmincon'
        ResultsVec = fmincon(SolverStruct.prData);
    case 'IPOPT'
        ! ./BLOM_Ipopt/BLOM_NLP
        ResultsVec = load('result.dat');
end


SolverResult = BLOM_ConvertVectorToStruct(ModelSpec.all_names,ResultsVec);