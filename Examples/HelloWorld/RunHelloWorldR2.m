clear

%% Convert model to optimization problem
[ModelSpec,block2,stepVars2,allVars2] = BLOM_ExtractModel2('HelloWorldNew',10,1,'Euler',1);
%%
[RunResults ResultsVec] = BLOM_RunModel(ModelSpec);
%%
SolverStruct = BLOM_ExportToSolver(ModelSpec,'IPOPT');

[OptGuess ExtVars InitialStates ] = BLOM_SplitResults(ModelSpec,RunResults);

SolverStructData =  BLOM_SetProblemData(SolverStruct,ModelSpec,OptGuess, ExtVars, InitialStates);

SolverResult  =  BLOM_RunSolver(SolverStructData,ModelSpec);


%% plot the results
subplot(211);
plot(SolverResult.System_state);

ylabel('X');

subplot(212);
plot(SolverResult.System_u);
ylabel('U');

