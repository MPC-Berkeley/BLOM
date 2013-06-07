

% reset the outout file
f = fopen('UT.out','wt');
fclose(f);

% Simple continous time model
discr = {'Euler','Trapez','RK4'};
solver = {'IPOPT','fmincon'}
for i_dis = 1:length(discr)
    for i_sol = 1:length(solver)
        
        ModelSpec = BLOM_ExtractModel('HelloWorldCont',10,1,discr{i_dis});
        
        RunResults = BLOM_RunModel(ModelSpec);
        
        SolverStruct = BLOM_ExportToSolver(ModelSpec,solver{i_sol});
        
        [OptGuess ExtVars InitialStates ] = BLOM_SplitResults(ModelSpec,RunResults);
        
        
        SolverStructData =  BLOM_SetProblemData(SolverStruct,ModelSpec,OptGuess, ExtVars, InitialStates);
        
        SolverResult  =  BLOM_RunSolver(SolverStructData,ModelSpec);
        
        f = fopen('UT.out','at');
        fprintf(f,'\nHelloWorldCont %s %s\n',discr{i_dis},solver{i_sol});
        fprintf(f,'%f ',SolverResult.System_Continuous_state);
        fclose(f);
    end
end

% Simple discrete model
solver = {'IPOPT','fmincon'};
for i_sol = 1:length(solver)
    
    ModelSpec = BLOM_ExtractModel('HelloWorld',10);
    
    RunResults = BLOM_RunModel(ModelSpec);
    
    SolverStruct = BLOM_ExportToSolver(ModelSpec,solver{i_sol});
    
    [OptGuess ExtVars InitialStates ] = BLOM_SplitResults(ModelSpec,RunResults);
    
    
    SolverStructData =  BLOM_SetProblemData(SolverStruct,ModelSpec,OptGuess, ExtVars, InitialStates);
    
    SolverResult  =  BLOM_RunSolver(SolverStructData,ModelSpec);
    
    f = fopen('UT.out','at');
    fprintf(f,'\nHelloWorld %s\n',solver{i_sol});
    fprintf(f,'%f ',SolverResult.System_xk1);
    fclose(f);
end

% larger model with sin/cos
y0=-20;
for i_sol = 1:length(solver)
    
    ModelSpec = BLOM_ExtractModel('SlidesExample',60,0.5);
    
    RunResults = BLOM_RunModel(ModelSpec);
    
    SolverStruct = BLOM_ExportToSolver(ModelSpec,solver{i_sol});
    
    [OptGuess ExtVars InitialStates ] = BLOM_SplitResults(ModelSpec,RunResults);
    
    SolverStructData =  BLOM_SetProblemData(SolverStruct,ModelSpec,OptGuess, ExtVars, InitialStates);
    
    SolverResult  =  BLOM_RunSolver(SolverStructData,ModelSpec);
    f = fopen('UT.out','at');
    fprintf(f,'\nSlidesExample %s \n',solver{i_sol});
    fprintf(f,'%f ',SolverResult.System_yk);
    fclose(f);
end

%%
clear 
clc
load RegularizationData

ModelSpec = BLOM_ExtractModel('Regularization',pred_horizon,dt); 

RunResults = BLOM_RunModel(ModelSpec);

SolverStruct = BLOM_ExportToSolver(ModelSpec,'IPOPT');

[OptGuess ExtVars InitialStates ] = BLOM_SplitResults(ModelSpec,RunResults);

SolverStructData =  BLOM_SetProblemData(SolverStruct,ModelSpec,OptGuess, ExtVars, InitialStates);

SolverResult  =  BLOM_RunSolver(SolverStructData,ModelSpec);
f = fopen('UT.out','at');
fprintf(f,'\nEPMO %s \n','IPOPT');
fprintf(f,'%f ',SolverResult.BilinearModel);
fclose(f);


