% Quick script to run SimpleNonLinProb.mdl

BLOM_SetDataLogging('SimpleNonLinProb')

ModelSpec = BLOM_ExtractModel('SimpleNonLinProb',1);

RunResults = BLOM_RunModel(ModelSpec);

SolverStruct = BLOM_ExportToSolver(ModelSpec,'IPOPT');

[OptGuess ExtVars InitialStates ] = BLOM_SplitResults(ModelSpec,RunResults);

SolverStructData =  BLOM_SetProblemData(SolverStruct,ModelSpec,OptGuess, ExtVars, InitialStates);

SolverResult  =  BLOM_RunSolver(SolverStructData,ModelSpec);