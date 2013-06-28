clear
clc
%% define variables
load model
T0=21; %current temp
mdotmax=740;
mdotmin=220;

simTime=12; %in hours
predHorizon=6; %in hours

weatherPred=smap_read_array('395005af-a42c-587f-9c46-860f3061ef0d',48.5,48.5);
weatherPred=(weatherPred-32)*5/9;
weatherPred=resample(weatherPred,weatherPred.Time(1):datenum(0,0,0,0,15,0):weatherPred.Time(end));
Time=weatherPred.Time(1:4*simTime);
weatherPred=arx_prep(weatherPred.Data,3);

%% Setup optimization problem
System_Load = 0;
for i=1:length(model)
    ModelSpec{i} = BLOM_ExtractModel('CITRIS2',1+4*predHorizon);
    RunResults{i} = BLOM_RunModel(ModelSpec{i});
    SolverStruct{i} = BLOM_ExportToSolver(ModelSpec{i},'IPOPT');
    [OptGuess{i} ExtVars InitialStates{i}] = BLOM_SplitResults(ModelSpec{i},RunResults{i});
end

for i=1:4*simTime %4*number of hours (model is discretized in 15 min intervals)
    %% run MPC
    ExtVars.System_Load=weatherPred(:,i:(i+4*predHorizon))';
        for j=1:length(model)
            SolverStructData{j} =  BLOM_SetProblemData(SolverStruct{j},ModelSpec{j},OptGuess{j}, ExtVars, InitialStates{j});
            try
                SolverResult{j}  =  BLOM_RunSolver(SolverStructData{j},ModelSpec{j});
            catch
                fprintf('I get here')
            end
            controlInputs{j}=[SolverResult{j}.System_mdot(1);SolverResult{j}.System_vlv(1)];
            controlInputHistory{j}(:,i)=controlInputs{j};
            tempHistory{j}(i)=InitialStates{j}.System_T;
            %% update currentTemp
            InitialStates{j}.System_T=model{j}.A*InitialStates{j}.System_T+model{j}.B1*controlInputs{j}(1)*controlInputs{j}(2)+model{j}.B2*controlInputs{j}(1)*InitialStates{j}.System_T+model{j}.W*weatherPred(:,i)+model{j}.C;
            OptGuess{j}=SolverResult{j};
        end
    disp('15 minutes has passed')
    disp('')
end

save simResults controlInputHistory tempHistory
%% plotting
mdot=controlInputHistory{1}(1,:)';
vlv=controlInputHistory{1}(2,:)';
temp=tempHistory{1}';

for i=2:length(model)
    mdot=[mdot controlInputHistory{i}(1,:)'];
    vlv=[vlv controlInputHistory{i}(2,:)'];
    temp=[temp tempHistory{i}'];
end

figure(4)
plot(Time,mdot)
datetick('x',15)
axis([Time(1) Time(end) mdotmin-10,mdotmax+10])
xlabel('Time')
ylabel('CFM')
figure(5)
plot(Time,vlv)
datetick('x',15)
axis([Time(1) Time(end) -10,110])
xlabel('Time')
ylabel('Valve position (%)')
figure(6)
plot(Time,temp,Time,22*ones(length(Time)),'-.',Time,20*ones(length(Time)),'-.')
datetick('x',15)
axis([Time(1) Time(end) 20-3,22+3])
xlabel('Time')
ylabel('Temperature (^oC)')