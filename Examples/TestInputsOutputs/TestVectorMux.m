
ModelSpec = BLOM_ExtractModel('VectorMux',1,[],[],[]);

RunResults = BLOM_RunModel(ModelSpec);

SolverStruct = BLOM_ExportToSolver(ModelSpec,'ipopt');

[OptGuess ExtVars InitialStates] = BLOM_SplitResults(ModelSpec,RunResults);

SolverStructData = BLOM_SetProblemData(SolverStruct, ModelSpec, OptGuess, ...
    ExtVars, InitialStates);

SolverResult = BLOM_RunSolver(SolverStructData,ModelSpec);
