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
        
        % mfilename('fullpath') returns the entire path to this script
        % the first output of fileparts(ans) gives the path string,
        % same as dirname in unix
        BLOM_dir = fileparts(mfilename('fullpath'));
        
        BLOM_NLP_exe = [ BLOM_dir '/BLOM_Ipopt/BLOM_NLP' ];
        if (~exist(BLOM_NLP_exe,'file'))
            error(['BLOM_NLP not found at ' BLOM_NLP_exe '. Run BLOM_Setup.']);
            SolverResult = [];
            ResultsVec   = [];
            return;
        end
        
        res = system(BLOM_NLP_exe);
        ResultsVec = load('result.dat');
end


SolverResult = BLOM_ConvertVectorToStruct(ModelSpec.all_names_struct,ResultsVec);