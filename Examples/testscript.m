% test script
simin = [0 1];

profile on;
[a,block,allVars] = BLOM_ExtractModel('testBFS',1,1,1,1);
profile viewer;

%% note, should find the following blocks in the BFS search (parentheses
%% has the number of outports that each block should have if it has more
%% than one)

testBFS_Blocks = {'Add1','Integrator','Add','InputFromWorkspace','InputFromSimulink',...
    'Polyblock2','ExternalFromWorkspace','ExternalFromSimulink1','Add3','Add2','Constant2','Constant3',...
    'Unit Delay1','Polyblock5','From','Unit Delay','Polyblock','From1','Demux','Product',...
    'Polyblock3','Subsystem','Subsystem/Unit Delay2','Subsystem/subadd1',...
    'Subsystem/Product1','Subsystem/Constant8','Subsystem/Constant5','From2','From3',...
    'Subsystem1','Subsystem1/sub_Polyblock1',...
    'Subsystem1/sub_Product1','Subsystem1/sub_Constant4','Subsystem1/sub_Constant5','Constant8'};

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