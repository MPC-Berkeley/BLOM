function [ SolverResult, ResultsVec, ResultInfo ] = ...
    BLOM_RunSolver(SolverStruct,ModelSpec,options)
%
%   [ SolverResult  ResultsVec ]=
%   BLOM_RunSolver(SolverStruct,ModelSpec,options)
%
%   Runs the solver and returns the solution .
%
% Input:
%   SolverStruct - Solver description struct, created with BLOM_ExportToSolver.
%   ModelSpec -    Model structure generatated by BLOM_ExtractModel.
%   options   -    options structure created by BLOM_optset or optimset function.
%
% Output:
%   SolverResult -  Structure with fields according to ModelSpec, holding
%                   the solver results. 
%   ResultsVec -    Vector with the same results    


switch lower(SolverStruct.solver)
    case 'fmincon'
        if nargin > 2
            SolverStruct.prData.options = optimset( ...
                SolverStruct.prData.options, options);
        end
        [ ResultsVec, FVAL,EXITFLAG] = fmincon(SolverStruct.prData);
    case 'linprog'
        if nargin > 2
            SolverStruct.prData.options = optimset( ...
                SolverStruct.prData.options, options);
        end
        [ ResultsVec, FVAL,EXITFLAG] = linprog(SolverStruct.prData);
    case 'ipopt'
        
        % mfilename('fullpath') returns the entire path to this script
        % the first output of fileparts(ans) gives the path string,
        % same as dirname in unix
        BLOM_dir = fileparts(mfilename('fullpath'));
        
        BLOM_NLP_exe = fullfile(BLOM_dir,'BLOM_Ipopt','BLOM_NLP');
        if ispc
            BLOM_NLP_exe = [BLOM_NLP_exe '.exe'];
        end
        if (~exist(BLOM_NLP_exe,'file'))
            error(['BLOM_NLP not found at ' BLOM_NLP_exe '. Run BLOM_Setup.']);
        end
        result_filename = fullfile(pwd,'result.dat');
        if exist(result_filename,'file')
            delete(result_filename) % delete any previous results file
        end
        EXITFLAG = system(BLOM_NLP_exe);
        if exist(result_filename,'file')
            ResultsVec = load(result_filename);
        else
            error('%s did not execute successfully, \n exit flag was %d', ...
                BLOM_NLP_exe, EXITFLAG)
        end
end


SolverResult = BLOM_ConvertVectorToStruct(ModelSpec.all_names_struct,ResultsVec);
ResultInfo.SolverExitFlag = EXITFLAG;
