% test script
simin = [0 1];

profile on;
a = BLOM_ExtractModel('testBFS',1,1,1,1);
profile viewer;

%% note, should find the following blocks in the BFS search (parentheses
%% has the number of outports that each block should have if it has more than one)

% Add1, Integrator, Add, InputFromWorkSpace, InputFromSimulink, Polyblock2 (2),
% ExternalFromWorkSpace, ExternalFromSimulink, Mux?, Add3, Add2, Constant2, Constant3,
% Unit Delay1, Polyblock5, Polyblock3, From, UnitDelay, Polyblock(2),
% Subsystem, (Inside Subsystem: Unit Delay2, subadd1, product1, Constant8,
% Constant5), Subsystem1 (2), (Inside Subsystem1: sub_Polyblock1, sub_Product1,
% sub_Constant4, sub_Constant5), Constant8
% TOTAL: 31 blocks to find