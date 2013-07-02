% test script
simin = [0 1];

%profile on;
horizon = 3;
[a,block,stepVars,allVars] = BLOM_ExtractModel('testBFS',horizon);
%profile viewer;

%% note, should find the following blocks in the BFS search (parentheses
%% has the number of outports that each block should have if it has more
%% than one)

testBFS_Blocks = {'Add1','Unit Delay2','Add','InputFromWorkspace','InputFromSimulink',...
    'Polyblock2','ExternalFromWorkspace','ExternalFromSimulink1','Add3','Add2','Constant3',...
    'Unit Delay1','Polyblock5','From','Unit Delay','Polyblock','From1','Demux','Product',...
    'Polyblock3','Subsystem','Subsystem/Unit Delay2','Subsystem/subadd1',...
    'Subsystem/Product1','Subsystem/Constant8','Subsystem/Constant5','Subsystem/In1','Subsystem/In2',...
    'From2','From3','From4','From5','Mux',...
    'Subsystem1','Subsystem1/sub_Polyblock1','Subsystem1/In1',...
    'Subsystem1/sub_Product1','Subsystem1/sub_Constant4','Subsystem1/sub_Constant5','Constant8',...
    'Subsystem2','Subsystem2/In1','Subsystem2/In2','Subsystem2/Add4','Subsystem2/Gain','Subsystem2/Integrator1',...
    'Subsystem2/Subsubsystem1', 'Subsystem2/Subsubsystem1/In1', 'Subsystem2/Subsubsystem1/Gain',...
    'Subsystem2/Subsubsystem2', 'Subsystem2/Subsubsystem2/In1', 'Subsystem2/Subsubsystem2/Gain',...
    'Bound', 'Bound1', 'DiscreteCost','Bound2','Demux1','Bound3'};

for i = 1:length(testBFS_Blocks)
    testBFS_Blocks{i} = ['testBFS/' testBFS_Blocks{i}];
end
    
% for testBFS, make sure that the proper blocks were found from the BFS

fprintf('-----------------------------------\n')
fprintf('-----------------------------------\n')

fprintf('Blocks incorrectly found:\n')
for i = 1:length(block.names)
    % test to see if it found incorrect blocks
    if any(strcmp(testBFS_Blocks,block.names{i}))
        % do nothing, it found the block
    else
        block.names{i}
    end
end

fprintf('\n\n-----------------------------------\n')
fprintf('-----------------------------------\n')
fprintf('Blocks that were not found:\n')
for i = 1:length(testBFS_Blocks)
    % test to see if it found all the correct blocks
    if any(strcmp(testBFS_Blocks{i},block.names))
        % do nothing, it found the block
    else
        testBFS_Blocks{i}
    end
end
