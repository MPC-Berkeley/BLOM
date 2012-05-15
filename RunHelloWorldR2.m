

clear

%% Convert model to optimization problem
ModelSpec = BLOM_ExtractModel('HelloWorldCont',10,1,'Euler');

RunResults = BLOM_RunModel(ModelSpec);


SolverStruct = BLOM_ExportToSolver(ModelSpec,'fmincon');

[OptGuess ExtVars InitialStates ] = BLOM_SplitResults(ModelSpec,RunResults);


SolverStructData =  BLOM_SetProblemData(SolverStruct,ModelSpec,OptGuess, ExtVars, InitialStates);

SolverResult  =  BLOM_RunSolver(SolverStructData,ModelSpec);


%% plot the results
subplot(211);
plot(SolverResult.System_Continuous_state);

ylabel('X');

subplot(212);
plot(SolverResult.System_u);
ylabel('U');

