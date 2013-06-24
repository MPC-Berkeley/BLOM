%> @file BLOM_ExtractModel.m
%> @brief BLOM 2.0 Outline based on what we discussed previously

%> Overview of BLOM Extract Model Structure
%> - Find all Input, External, Constraint and Cost BLOM library blocks.
%> - Build cost and constraint dependency wire graph, starting from Constraint and Cost blocks and traverse backwards.
%>     -  BFS to find the graph.
%>     -  Recognize reaching non-BLOM inputs (user error).
%>     -  Format of the graph: Sparse connectivity matrix of all wires, 
%>         store just the block handles.
%>     -  Each wire has a property whether optimization variables is
%>             introduced for it for each derivative evaluation, 
%>         intermediate major time step, initial time, final time.
%>     -  Store for each wire, if it is constrained, cost (discrete or continuous), 
%>         input, external.
%> - Flatten the wires:
%>     -  "Fetch" all PolyBlocks data: variables from the base workspace or values from the dialog (just do evalin('base',get_param(h,'A')).
%>     -  Get all possible dimensions (with or without the Simulink "compile"). 
%>     -  Translate all supported Simulink blocks to PolyBlocks using the known dimensions.   
%>     -  Eliminate duplicate routing variables. 
%>     -  Combine all PolyBlocks to a large P and K matrices.
%>     -  Eliminate whatever is possible (constants, identity PolyBlocks, linearly dependent constraints, replicated P rows, unused terms).
%>     -  Propagate time property information using the processed P and K structure.    
%> - Create the whole problem:
%>     -  Introduce and duplicate all major time step variables, including their dependencies.
%>     -  Introduce and duplicate minor time steps variables and their dependencies.
%>     -  Link by the appropriate Butcher tableau values major and minor time steps variables.
%>     -  Add move blocking constraints.
%>     -  Eliminate whatever is possible.
%>     -  Mark external, input.
%>     -  Take care of cost.

%======================================================================
%> @brief This function extracts a Simulink model and is the first step
%> into converting it into an optimization problem
%>
%> Examples:
%>   ModelSpec = BLOM_ExportModel; 
%>       Converts current active model.
%>
%>   ModelSpec = BLOM_ExportModel('Sys',10,0.5);  
%>       Converts current model 'Sys', with a a time step of 0.5 sec and 10 time steps.
%>   ModelSpec = BLOM_ExportModel('Sys',10,0.5,);  
%>       Converts current model 'Sys', with a a time step of 0.5 sec and 10 time steps, 
%>       assuming discrete model.
%>
%>   ModelSpec = BLOM_ExportModel('Sys',10,0.5,'Trapez');  
%>       Converts current model 'Sys', with a a time step of 0.5 sec ,10
%>       time steps and trapezoidal discretization method.
%>
%> @param name name of the model to convert
%> @param horizon integer number of steps in prediction horizon length, 1 if not
%> supplied
%> @param dt time step size [s], 1 if not supplied
%> @param disc discretization method {'none,'Euler','Trapez','RK4'}
%> @param options options created by BLOM_opset function
%>
%> @retval out1 return value for the first output variable
%> @retval out2 return value for the second output variable
%======================================================================

function [ModelSpec,block,stepVars,allVars] = BLOM_ExtractModel2(name,horizon,dt,integ_method,options)
    % load system. does nothing if model is not open
    load_system(name);
    if ~exist('dt','var') || isempty(dt)
        dt = 1;
    end
    if ~exist('integ_method','var') || isempty(integ_method)
        integ_method = 'None';
    end
    if ~exist('options','var')
        options = [];
    end
    % evaluate model to get dimensions
    eval([name '([],[],[],''compile'');']); 
    try
        [boundHandles,costHandles,inputAndExternalHandles] = findBlocks(name);
        [block,stepVars,stop] = ...
            searchSources(boundHandles,costHandles,inputAndExternalHandles,name);
        stepVars = setBounds(block, stepVars);
        stepVars = setCosts(block, stepVars);
        % FIX: should implement something that stops the code after analyzing
        % all the blocks and finding an error in the structure of the model
        
        % expands all the inport and outport indexes such that it points to
        % each scalar instead of just to first entry
        block = expandInportOutportIdx(block);
        
        if stop == 1
            % break the code somehow?
        end
        
        stepVars = labelTimeRelevance(stepVars,block,inputAndExternalHandles);
                
        % create allVars. allVars contains all variables at all time steps
        % and is made from stepVars
        allVars = createAllVars(stepVars,horizon);
        
        block = expandBlock(block,horizon,stepVars,allVars);

        
        stepVars.optVarIdx = cleanupOptVarIdx(stepVars.optVarIdx);
        
        allVars = allOptVarIdxs(allVars,block,stepVars,horizon);
        
        % create large P and K matrix for entire problem
        % initially we create the P and K matrix of stepVars then we create
        % the P and K matrix for allVars
        try
           [stepP,stepK] = combinePK(block,stepVars);
        catch err
           rethrow(err)
        end
        
        % create P and K matrix for allVars
        [allP,allK] = createAllPK(stepP,stepK,stepVars,horizon,allVars);
        
        % convert to ModelSpec - this part will merge BLOM 2.0 and BLOM 1.0
        ModelSpec = convert2ModelSpec(name,horizon,integ_method,dt,options,stepVars,allVars,block, allP,allK);
        % following code checks whether or not inports and outportHandles
        % was filled in properly
        fprintf('--------------------------------------------------------\n')
        fprintf('Test block.stepOutputIdx and inports is correct\n')
        fprintf('Everything is good if nothing prints\n')
        fprintf('--------------------------------------------------------\n')
        for i = 1:length(block.handles)
            currentPorts = get_param(block.handles(i),'PortHandles');
            currentOutports = currentPorts.Outport;
            currentInports = currentPorts.Inport;
            % compare lists of Outports
            
            currentBlockOutports = zeros(length(block.stepOutputIdx{i}),1);
            for idx = 1:length(block.stepOutputIdx{i})
                currentBlockOutports(idx) = stepVars.outportHandle(block.stepOutputIdx{i}(idx));
            end
            
            diffOutports = setdiff(currentBlockOutports,currentOutports);
            if ~isempty(diffOutports)
                fprintf('Difference in Outports in %s\n',block.names{i})
                currentOutports
                get_param(stepVars.outportHandle(block.stepOutputIdx{i}),'Parent')
                fprintf('--------------------------------------------------------\n')
            end
            
            % compare inports
            inportsOutport = zeros(length(currentInports),1);
            for index = 1:length(currentInports)
                currentLine = get_param(currentInports(index),'Line');
                if currentLine > 0
                    inportsOutport(index) = get_param(currentLine,'SrcPortHandle');
                end
            end
            
            currentBlockInports = zeros(length(block.stepInputIdx{i}),1);
            for idx = 1:length(block.stepInputIdx{i})
                currentBlockInports(idx) = stepVars.outportHandle(block.stepInputIdx{i}(idx));
            end
            
            diffInports = setdiff(inportsOutport,currentBlockInports);
            if ~isempty(diffInports) && ~any(block.handles(i)==inputAndExternalHandles)
                fprintf('Difference in inputs in %s\n',block.names{i})
                ['inputsOutport;' inportsOutport]
                block.stepInputIdx{i}
                get_param(inportsOutport,'Parent')
                fprintf('--------------------------------------------------------\n')
            end
            
            
            % first need to get list of inports for each block and the
            % outport that's associated with that inport
            
        end
        
        %following code is to make sure searchSources works
        fprintf('The Number of blocks found is %.0f\n',length(block.handles))
        fprintf('The Number of outports found is %.0f\n',length(stepVars.outportHandle))
        for i = 1:length(stepVars.outportHandle);
            parent = get_param(stepVars.outportHandle(i),'Parent');
            stepVars.outportHandle(i);
            portType = get_param(stepVars.outportHandle(i),'PortType');
            if ~strcmp(portType,'outport')
                fprintf('Oops, this is not an outport, look to see what happened\n')
                portType
            end
        end
        

        fprintf('\n\n\n Here we have the stored data from block\n\n\n');

        for i = 1:length(block.handles);
            someBlock = block.names{i};
        end
        
        % check to see if stepVars points to the proper block
        fprintf('Test to make sure stepVars points to the proper block\n')
        fprintf('--------------------------------------------------------\n')
        for i = 1:length(stepVars.block)
            currentBlockName = block.names(stepVars.block(i));
            parent = get_param(stepVars.outportHandle(i),'Parent');
            
            if strcmp(currentBlockName,parent)
                %fprintf('correct index\n')
            else
                fprintf('INCORRECT INDEX, BAD\n')
                currentBlockName
                parent
                fprintf('\nExambine Above\n')
            end
        end
        fprintf('--------------------------------------------------------\n')

    catch err
        % close evaluation of models
        eval([name '([],[],[],''term'');']);
        rethrow(err)
    end
    % close evaluation of models
    eval([name '([],[],[],''term'');']);
    
end

%%
%======================================================================
%> @brief Finds all Input, External, Constraint and Cost BLOM library
%> blocks
%>
%> More detailed description of the problem.
%>
%> @param name Name of the simulink model
%>
%> @retval blockHandles handles of all the relevant blocks
%======================================================================

function [boundHandles,costHandles,inputAndExternalHandles] = findBlocks(name)
    % FindAll may not be the most efficient way to find the handles
    boundHandles = [find_system(name,'FindAll','on','ReferenceBlock','BLOM_Lib/Bound')];
    costHandles = [find_system(name,'FindAll','on','ReferenceBlock','BLOM_Lib/DiscreteCost')];
    inputAndExternalHandles = ...
        [find_system(name,'FindAll','on','ReferenceBlock','BLOM_Lib/InputFromWorkspace'); ...
        find_system(name,'FindAll','on','ReferenceBlock','BLOM_Lib/InputFromSimulink'); ...
        find_system(name,'FindAll','on','ReferenceBlock','BLOM_Lib/ExternalFromWorkspace'); ...
        find_system(name,'FindAll','on','ReferenceBlock','BLOM_Lib/ExternalFromSimulink')];
end

%%
%======================================================================
%> @brief Find all connected source blocks using a breadth first search.
%> Includes the handles already given.
%>
%> Also takes into consideration what type of block it's looking at
%>
%> @param boundHandles array of handles of all the bound blocks. This is where we
%> begin the BFS
%> @param costHandles array of handles all of the cost blocks. We also
%> begin the BFS from these blocks.
%> @param varargin the external and input from simulink handles. the
%> sources of these blocks are not relevant to the optimization problem.
%> These are used to find out where a BFS branch should stop.
%>
%> @retval block structure with information about all the blocks found in
%> the BFS
%> @retval stepVars structure that contains information about each
%> optimization variable
%> @retval stop if there are any unconnected blocks or blocks that BLOM
%> does not support, this parameter gives a 1 to indicate true and 0 to
%> indicate false
%======================================================================

function [block,stepVars,stop] = searchSources(boundHandles,costHandles,...
    fromSimulink,name)
    % only flag stop = 1 when there is a bad block
    stop = 0;
    
    % gets blocks of final stopping places
    fromSimulinkPorts = get_param(fromSimulink,'PortHandles');
    sOutportHandles = zeros(length(fromSimulinkPorts),1);
    sOutportIdx = 1;
    % add all the outport handles of the fromSimulink blocks. The BFS
    % will stop at these handles
    for idx = 1:length(fromSimulinkPorts);
        if length(fromSimulinkPorts)==1
            fromSimulinkPorts = {fromSimulinkPorts};
        end
        currentOutports = fromSimulinkPorts{idx}.Outport;
        sOutportLength = length(currentOutports);
        sOutportHandles(sOutportIdx:(sOutportIdx+sOutportLength-1)) = ...
            currentOutports;
        sOutportIdx = sOutportIdx+sOutportLength;
    end
    
    
    initialSize = 20;
    
    % initialize structures and fields
    
    % stepVars stores information for every optimization variable
    stepVars.zeroIdx = 1; % index of first stepVars zero
    stepVars.block = zeros(initialSize,1); % the index of the block
    stepVars.outportNum = zeros(initialSize,1); % outport number. e.g. if there are multiple outports of a block
    stepVars.outportHandle = zeros(initialSize,1); % handle of specific outport
    stepVars.outportIndex = zeros(initialSize,1); % index of specific outports. normally just 1 (will be more than 1 if the outport is a vector)
    stepVars.optVarIdx = (1:initialSize)'; % points to the true optimization variable index. (will usually point to itself, otherwise, some "true" variable)
    stepVars.initCost = zeros(initialSize,1); %default cost is zero. 1 if we see that it's part of the cost at initial time step
    stepVars.interCost = zeros(initialSize,1); %default cost is zero. 1 if we see that it's part of the cost at intermediate time steps
    stepVars.finalCost = zeros(initialSize,1); %default cost is zero. 1 if we see that it's part of the costat final time step
    stepVars.initUpperBound = inf*ones(initialSize,1); % upper bound of at initial time step. default inf
    stepVars.initLowerBound = -inf*ones(initialSize,1); % lower bound of at initial time step. default -inf
    stepVars.interUpperBound = inf*ones(initialSize,1); % upper bound of at intermediate time step. default inf
    stepVars.interLowerBound = -inf*ones(initialSize,1); % lower bound of at intermediate time step. default -inf
    stepVars.finalUpperBound = inf*ones(initialSize,1); % upper bound of at final time step. default inf
    stepVars.finalLowerBound = -inf*ones(initialSize,1); % lower bound of at final time step. default -inf
    stepVars.state = false(initialSize,1); %indicator as to whether the variable is a state variable
    stepVars.input = false(initialSize,1); %indicator as to whether the variable is an input variable
    stepVars.external = false(initialSize,1); %indicator as to whether the variable is an external variable
    % the following stepVars fields state whether or not a variable is
    % relevant at that time step. They are not filled in but put here for
    % reference
    %stepVars.initTime
    %stepVars.interTime
    %stepVars.finalTime
    
    % block stores information about each block
    block.zeroIdx = 1; % index of first block zero
    block.names = cell(initialSize,1); % name of block
    block.handles = zeros(initialSize,1); % handle of block
    block.P = cell(initialSize,1); % P matrix of block, if relevant
    block.K = cell(initialSize,1); % K matrix of block, if relevant
    block.stepInputIdx = cell(initialSize,1); % the stepVar indices of the outports connected to the inports of that block
    block.stepOutputIdx = cell(initialSize,1); % the stepVar indices of the outports of that block
    block.dimensions = cell(initialSize,1); % dimensions of each outports. first value is outport #, then second two values are dimensions of outport
    block.bound = false(initialSize,1); % indicator to whether block is a bound block
    block.cost = false(initialSize,1); % indicator to whether block is a cost block
    block.subsystem = false(initialSize,1);% indicator to whether block is a subsystem block and not a BLOM block
    block.mux = false(initialSize,1);% indicator to whether block is a mux block
    block.demux = false(initialSize,1);% indicator to whether block is a demux block
    block.delay = false(initialSize,1);% indicator to whether block is a delay block
    block.subsystemInput = false(initialSize,1);% indicator to whether block is a input block within subsystem
    block.inputBlock = false(initialSize,1); % indicator for whether this is a InputFromSimulink OR InputFromWorkspace
    block.externalBlock = false(initialSize,1); % indicator for whether this is an ExternalFromWorkSpace OR ExternalFromSimulink
    block.fromBlock = false(initialSize,1);% indicator to whether block is a from block
    block.blockType = cell(initialSize,1);% store the blockType for block for BFS switching
    block.refBlock = cell(initialSize,1);% store the refblock for block for BFS switching
    block.reroute = false(initialSize,1);% true if this block will need to reroute variables in one way or another
    block.specialFunPresence = false(initialSize,1);% true if this block contains a special function


    
    % find all lines connected to costs and bounds and then get outport
    % ports from there. fill in stepVars and block structures as necessary
    
    % iZero is index of first zero of outportHandles
    outports.outportHandles = zeros(initialSize,1);
    outports.zeroIdx = 1; % the index of the first zero of this outportHandles
    outports.searchIdx = 1; % the current index we're looking at in the BFS
    
    % get all outports connected to the bounds 
    for i = 1:length(boundHandles)
        portH = get_param(boundHandles(i),'PortHandles');
        % costs and bounds should only have one inport
        currentInport = [portH.Inport];
        [block, currentBlockIndex] = updateBlock(block,currentInport);
        [outports,~] = addToBFS(outports,currentInport,block,currentBlockIndex);
    end
    
    % get all outports connected to costs
    for i = 1:length(costHandles)
        portH = get_param(costHandles(i),'PortHandles');
        % costs and bounds should only have one inport
        currentInport = [portH.Inport];
        [block, currentBlockIndex] = updateBlock(block,currentInport);
        [outports,~] = addToBFS(outports,currentInport,block,currentBlockIndex);
    end
 
    while 1
        if outports.outportHandles(outports.searchIdx) == 0
            % all handles have been found and searched
            break
        elseif any(outports.outportHandles(outports.searchIdx)==sOutportHandles)
            % if the current handle is equal to any of the external or input
            % from simulink, then continue with the loop and do not look
            % for more branches from that block
            % in this case do we not need to run addToBFS
            currentOutport = outports.outportHandles(outports.searchIdx);
            [block, currentBlockIndex] = updateBlock(block,currentOutport);
            [stepVars,block] = updateStepVars(stepVars,block,currentBlockIndex,currentOutport);
            outports.searchIdx = outports.searchIdx + 1;
            continue
        elseif outports.searchIdx> length(outports.outportHandles)
            % should not reach here, put here just in case
            break
        elseif ~outports.outportHandles(end) == 0
            % no more space for more handles, allocate more space
            outports.outportHandles = [outports.outportHandles;...
                zeros(length(outports.outportHandles)*2,1)];
        end
        
        % call updateBlock here somehow
        % call addtobfs here
        % rerouting case
        currentOutport = outports.outportHandles(outports.searchIdx);
        [block, currentBlockIndex] = updateBlock(block,currentOutport);
        [outports,allOutportsFound] = addToBFS(outports,currentOutport,block,currentBlockIndex);
        
        if block.reroute(currentBlockIndex)
            % this block will need some rerouting. mux and demux will
            % require special rerouting
            if block.mux(currentBlockIndex)
                % special rerouting for mux
                % the original variables are the outports connected to the
                % mux
                mux_optVarIdx = zeros(length(allOutportsFound),1);
                for idx=1:length(allOutportsFound)
                    tempOutport = allOutportsFound(idx);
                    [block,tempBlockIndex] =...
                        updateBlock(block,tempOutport);
                    [stepVars,block,sameOptIndex] =...
                        updateStepVars(stepVars,block,tempBlockIndex,tempOutport);
                    mux_optVarIdx(idx) = sameOptIndex;
                end
                
                [stepVars,block] = updateStepVars(stepVars,block,currentBlockIndex,currentOutport,mux_optVarIdx);
                
            elseif block.demux(currentBlockIndex)
                % special rerouting for demux
                tempOutport = allOutportsFound;
                [block,tempBlockIndex] =...
                    updateBlock(block,tempOutport);
                [stepVars,block,sameOptIndex] =...
                    updateStepVars(stepVars,block,tempBlockIndex,tempOutport);
                
                demuxPorts = get_param(block.handles(currentBlockIndex),'PortHandles');
                demuxOutports = [demuxPorts.Outport];
                
                for ii = 1:length(demuxOutports)
                    [stepVars,block] = updateStepVars(stepVars,block,currentBlockIndex,currentOutport,sameOptIndex+ii-1);
                end
                
            else
                % rerouting procedure for Subsystems (going in and out),
                % from/goto blocks. There should only be one outport
                [block, tempBlockIndex] = updateBlock(block,allOutportsFound);
                % first fill in stepVars for the "original" variable
                [stepVars,block,sameOptIndex] = updateStepVars(stepVars,block,tempBlockIndex,allOutportsFound);
                % then fill in the rerouting 
                [stepVars,block] = updateStepVars(stepVars,block,currentBlockIndex,currentOutport,sameOptIndex);
            end
        else
            % all other blocks
            [stepVars,block] = updateStepVars(stepVars,block,currentBlockIndex,currentOutport);
        end
        
        outports.searchIdx = outports.searchIdx + 1;

    end
    
    
    
    if any(outports.outportHandles==0)
        outports.outportHandles = outports.outportHandles(1:outports.zeroIdx-1);
    end
    
    block = updateInputsField(block,stepVars,outports.outportHandles);
    
    %check if any handles are -1 and remove them
    if any(outports.outportHandles==-1)
        % FIX PRINT STATEMENT. NOT NECESSARILY TRUE
        fprintf('One or more of the blocks are missing a connection\n');
        stop = 1;
    end
    
    % remove empty and 0 entries in blocks
    for field={'names', 'P','K','stepInputIdx','stepOutputIdx','dimensions',...
            'handles', 'bound', 'cost', 'subsystem', 'mux', 'demux', 'delay',...
            'subsystemInput', 'fromBlock', 'inputBlock', 'externalBlock', ...
            'blockType', 'refBlock', 'reroute', 'specialFunPresence'}
        block.(field{1}) = block.(field{1})(1:(block.zeroIdx-1));
    end

    
    % remove empty and 0 entries in stepVars
    for field={'block', 'outportNum','outportHandle','outportIndex',...
        'optVarIdx','initCost','interCost','finalCost','state','input','external','initLowerBound','initUpperBound',...
        'interLowerBound','interUpperBound','finalLowerBound','finalUpperBound'}
        stepVars.(field{1}) = stepVars.(field{1})(1:(stepVars.zeroIdx-1));
    end

    % FIX: need to find some way to remove all -1 handles. using setdiff
    % with [-1] reorders all the outport handles and puts it in ascending
    % order
    % outportHandles = setdiff(outportHandles,[-1]);
    
end

%======================================================================
%> @brief From given inports, see which outports are relevant and update
%> the list of outports accordingly
%>
%> More detailed description of the problem.
%>
%> @param inports inports of current source
%> @param existingOutports current outports found. used to check for
%> redundancy
%> @param iZero current index of the first zero 
%>
%> @retval outportHandles new outports found
%> @retval iZero updates current index of the first zero
%======================================================================

function [outports,allOutportsFound] = addToBFS(outports,inports,block,currentBlockIndex)
    
    currentOutport = outports.outportHandles(outports.searchIdx);
    
    if block.subsystem(currentBlockIndex) %FIX: may want to look into special BLOM blocks here
        % if there's a subsystem, inports is actually an array of the
        % outports of subsystem
        parent = get_param(currentOutport,'Parent');
        index = get_param(currentOutport, 'PortNumber');
        outportBlocks = find_system(parent,'SearchDepth',1,'regexp','on','BlockType','Outport');
        handle = get_param(outportBlocks{index},'Handle');
        portH = get_param(handle,'PortHandles');
        inports = [portH.Inport];
    elseif block.fromBlock(currentBlockIndex)
        tag = get_param(block.handles(currentBlockIndex),'GoToTag');
        name = get_param(block.handles(currentBlockIndex),'Parent');
        % there should only be one goto block associated with this from
        % block
        gotoBlock = find_system(name,'BlockType','Goto','GotoTag',tag);
        gotoPorts = get_param(gotoBlock{1},'PortHandles');
        % note: there should only be one outport associated with each goto
        % block
        inports = [gotoPorts.Inport];
    elseif block.subsystemInput(currentBlockIndex)
        parentOfBlock = get_param(block.handles(currentBlockIndex),'Parent');
        sourcePorts = get_param(parentOfBlock,'PortHandles');
        portNumber = eval(get_param(block.handles(currentBlockIndex),'Port'));
        inports = [sourcePorts.Inport];
        inports = inports(portNumber);
    else
        % from outport, find relevant inports
        parent = get_param(inports ,'Parent');
        portHandles = get_param(parent, 'PortHandles');
        inports = [portHandles.Inport];        
    end
   
    % found outports connected to inports provided
    allOutportsFound=zeros(length(inports),1);
    for i = 1:length(inports);
        currentLine = get_param(inports(i),'Line');
        % this gives the all the outports connected to this line
        outportFound = get_param(currentLine,'SrcPorthandle');
        allOutportsFound(i) = outportFound;
        outLength = length(outportFound);
        % in case outportHandles is too short
        if outLength > length(outports.outportHandles)-outports.zeroIdx+1;
            if outLength > length(outports.outportHandles)
                outports.outportHandles = [outports.outportHandles; zeros(outLength*2,1)];
            else
                outports.outportHandles = [outports.outportHandles; ...
                    zeros(length(outports.outportHandles),1)];
            end
        end
        
        diffOutports = setdiff(outportFound,outports.outportHandles);
        diffLength = length(diffOutports);
        outports.outportHandles(outports.zeroIdx:(diffLength+outports.zeroIdx-1))...
            = diffOutports;
        outports.zeroIdx = outports.zeroIdx + diffLength;
    end

end


%%
%======================================================================
%> @brief update stepVars structure 
%>
%> Description of how we use sameOptVar 
%> (or which outport we list as the original variable):
%> 1. Inports (within subsystems): the original variable the outport of the block 
%> that goes into the inport
%> 2. Outports (within subsystems): the original variable is the outport of
%> the block within the subsystem 
%> 3. From, GoTo Blocks: The original variable is the outport that connects
%> to the proper GoTo block
%> For mux and demux blocks, we assume that there's
%> either a scalars -> vector conversion or a vector -> scalars conversion. 
%> 4. Demux: The original variable is the outport connected to the demux
%> and that is being "demux'ed "
%> 5. Mux: The original variables are the outports connected to the mux
%> 
%> 
%> @param stepVars stepVars structure
%> @param stepVars.zeroIdx current index of first zero of stepVars
%> @param block.zeroIdx current index of first zero of block structure. used to
%> figure out which block the outport is from
%> @param currentOutport outport that you want to add information for
%> @param state lets the function know how to populate stepVars. current
%> states are 'bound', 'cost'
%>
%> @retval stepVars stepVars structure
%> @retval stepVars.zeroIdx updated index of first zero of stepVars.zeroIdx
%======================================================================

function [stepVars,block,varargout] = updateStepVars(stepVars,...
    block,currentBlockIndex,currentOutport,varargin)
%this function populates stepVars structure
    dimension = get_param(currentOutport,'CompiledPortDimensions');
    lengthOut = dimension(1)*dimension(2);
    
    if sum(any(stepVars.outportHandle==currentOutport)) == 0
        %only add general information for outports which haven't been looked at
        currentIndices = 1:lengthOut;
        portNumber = get_param(currentOutport,'PortNumber');
        
        % update the size of stepVars as necessary
        newLength = stepVars.zeroIdx+lengthOut;
        oldLength = length(stepVars.block);
        if newLength >= oldLength
            if newLength >= 2*oldLength
            % new entries greater than twice the old length
                for field={'block', 'outportNum','outportHandle','outportIndex',...
                        'cost'}
                        stepVars.(field{1}) = [stepVars.(field{1}); zeros(newLength*2,1)];
                end
                
                for field={'state', 'input','external'}
                        stepVars.(field{1}) = [stepVars.(field{1}); false(newLength*2,1)];
                end
                stepVars.optVarIdx = [stepVars.optVarIdx; ((oldLength+1):(newLength*2+oldLength))'];
                
                stepVars.initLowerBound = [stepVars.initLowerBound; -inf*ones(newLength*2,1)];
                stepVars.initUpperBound = [stepVars.initUpperBound; inf*ones(newLength*2,1)];
                stepVars.interLowerBound = [stepVars.interLowerBound; -inf*ones(newLength*2,1)];
                stepVars.interUpperBound = [stepVars.interUpperBound; inf*ones(newLength*2,1)];
                stepVars.finalLowerBound = [stepVars.finalLowerBound; -inf*ones(newLength*2,1)];
                stepVars.finalUpperBound = [stepVars.finalUpperBound; inf*ones(newLength*2,1)];
                
                stepVars.time = [stepVars.time; cell(newLength*2,1)];
            else % double the length of all fields in stepVars
                for field={'block', 'outportNum','outportHandle','outportIndex',...
                        'initCost', 'interCost', 'finalCost'}
                        stepVars.(field{1}) = [stepVars.(field{1}); zeros(oldLength,1)];
                end
                stepVars.optVarIdx = [stepVars.optVarIdx; ((oldLength+1):(oldLength*2))'];
                stepVars.initLowerBound = [stepVars.initLowerBound; -inf*ones(oldLength*2,1)];
                stepVars.initUpperBound = [stepVars.initUpperBound; inf*ones(oldLength*2,1)];
                stepVars.interLowerBound = [stepVars.interLowerBound; -inf*ones(oldLength*2,1)];
                stepVars.interUpperBound = [stepVars.interUpperBound; inf*ones(oldLength*2,1)];
                stepVars.finalLowerBound = [stepVars.finalLowerBound; -inf*ones(oldLength*2,1)];
                stepVars.finalUpperBound = [stepVars.finalUpperBound; inf*ones(oldLength*2,1)];
                
                for field={'state', 'input', 'external'}
                        stepVars.(field{1}) = [stepVars.(field{1}); false(oldLength,1)];
                end

            end
        end
        stepVars.block(stepVars.zeroIdx:(stepVars.zeroIdx+lengthOut-1)) = currentBlockIndex;
        stepVars.outportNum(stepVars.zeroIdx:(stepVars.zeroIdx+lengthOut-1)) = portNumber;
        stepVars.outportHandle(stepVars.zeroIdx:(stepVars.zeroIdx+lengthOut-1)) = currentOutport;
        stepVars.outportIndex(stepVars.zeroIdx:(stepVars.zeroIdx+lengthOut-1)) = currentIndices;
        
        block.stepOutputIdx{currentBlockIndex}(portNumber) = stepVars.zeroIdx;
        varIdx = stepVars.zeroIdx;
        stepVars.zeroIdx = stepVars.zeroIdx+lengthOut;
    else
        varIdx = find(stepVars.outportHandle==currentOutport,1);
    end
    
    if nargout == 3
        % if the nargout = 3, that means that we need to output the
        % optVarIdx for rerouting
        varargout{1} = varIdx;
    end
    

    % Here we do the switching for different types of blocks
    
    if block.bound(currentBlockIndex)
        % should never be here
       warning('Adding variable to bound block...WHY?')
        
    elseif block.cost(currentBlockIndex)
        % fill in cost field
        costHandle = block.handles(currentBlockIndex);
        stepVars.initCost(varIdx:(varIdx+lengthOut-1)) = strcmp(get_param(costHandle, 'initial_step'), 'on');
        stepVars.interCost(varIdx:(varIdx+lengthOut-1)) = strcmp(get_param(costHandle, 'intermediate_step'), 'on');
        stepVars.finalCost(varIdx:(varIdx+lengthOut-1)) = strcmp(get_param(costHandle, 'final_step'), 'on');
    elseif block.inputBlock(currentBlockIndex)
        % fill in input field
        stepVars.input(varIdx:(varIdx+lengthOut-1)) = true;
    elseif block.externalBlock(currentBlockIndex)
        % fill in external field
        stepVars.external(varIdx:(varIdx+lengthOut-1)) = true;
    elseif block.fromBlock(currentBlockIndex) || block.subsystem(currentBlockIndex) || block.subsystemInput(currentBlockIndex)...
            || block.demux(currentBlockIndex)
        % for from, subsystem and subsystem inport blocks, "point" to the
        % original outport
        if ~isempty(varargin)
            sameOptIndex = varargin{1};
            stepVars.optVarIdx(varIdx:(varIdx+lengthOut-1)) = ...
                (sameOptIndex):(sameOptIndex+lengthOut-1);
        end
    elseif block.delay(currentBlockIndex)
        % Unit Delay (or state case)
        stepVars.state(varIdx:(varIdx+lengthOut-1)) = true;
    elseif block.mux(currentBlockIndex)
        % Mux rerouting
        if ~isempty(varargin)
            sameOptIndex = varargin{1};
            stepVars.optVarIdx(varIdx:(varIdx+lengthOut-1)) = sameOptIndex;
        end
    else 
        % do nothing for all non special cases
        
    end
    
        
end

%%
%======================================================================
%> @brief given block structure and outport, populate the relevant fields
%> of block using currentOutport information
%>
%> More detailed description of the problem.
%>
%> @param block block structure
%> @param block.zeroIdx current index of first zero of block
%> @param currentOutport outport that you want to add information for
%>
%> @retval block block structure
%> @retval block.zeroIdx updated index of first zero of block
%======================================================================

function [block,currentBlockIndex] = updateBlock(block,currentOutport)

    currentBlockHandle = get_param(currentOutport,'ParentHandle');
    sourcePorts = get_param(currentBlockHandle,'PortHandles');
    
    % if the size of block equals the block.zeroIdx, need to double to
    % length of block
    if block.zeroIdx == length(block.handles)
        for field={'names', 'P','K','stepInputIdx','stepOutputIdx','dimensions','blockType','refBlock'}
                block.(field{1}) = [block.(field{1}); cell(block.zeroIdx,1)];
        end
        block.handles = [block.handles; zeros(block.zeroIdx,1)];
        block.bound = [block.bound; false(block.zeroIdx,1)];
        block.cost = [block.cost; false(block.zeroIdx,1)];
        block.subsystem = [block.subsystem; false(block.zeroIdx,1)];
        block.mux = [block.mux; false(block.zeroIdx,1)];
        block.demux = [block.demux; false(block.zeroIdx,1)];
        block.delay = [block.delay; false(block.zeroIdx,1)];
        block.reroute = [block.reroute; false(block.zeroIdx,1)];
        block.subsystemInput = [block.subsystemInput; false(block.zeroIdx,1)];
        block.fromBlock = [block.fromBlock; false(block.zeroIdx,1)];
        block.inputBlock = [block.inputBlock; false(block.zeroIdx,1)];
        block.externalBlock = [block.externalBlock; false(block.zeroIdx,1)];
        block.specialFunPresence = [block.specialFunPresence; false(block.zeroIdx,1)];
    end

    if ~any(block.handles==currentBlockHandle)
        % no duplicate blocks, add this block
        
        currentBlockName = get_param(currentOutport,'Parent');
        block.names{block.zeroIdx} = currentBlockName;
        
        % add the block type and reference block for specific cases
        block.blockType{block.zeroIdx} = get_param(currentBlockHandle,'BlockType');
        block.refBlock{block.zeroIdx} = get_param(currentBlockHandle,'ReferenceBlock');
        block.handles(block.zeroIdx) = currentBlockHandle;
        if strcmp(block.refBlock{block.zeroIdx},'BLOM_Lib/Polyblock')
            % store P&K matricies if the current block is a polyblock
            block.P{block.zeroIdx} = eval(get_param(currentBlockHandle,'P'));
            block.K{block.zeroIdx} = eval(get_param(currentBlockHandle,'K'));
            
            % we have to modify the P and K matrices because the columns of
            % P are only the inputs of the block so far. we need to include
            % the outputs
            totalOutputs = size(block.K{block.zeroIdx},1);
            block.P{block.zeroIdx} = blkdiag(block.P{block.zeroIdx},speye(totalOutputs));
            block.K{block.zeroIdx} = horzcat(block.K{block.zeroIdx},-speye(totalOutputs));
        elseif strcmp(block.refBlock{block.zeroIdx}, 'BLOM_Lib/Bound')
            block.bound(block.zeroIdx) = true;
        elseif strcmp(block.refBlock{block.zeroIdx}, 'BLOM_Lib/DiscreteCost')
            block.cost(block.zeroIdx) = true;
        elseif isempty(block.refBlock{block.zeroIdx}) && strcmp(block.blockType{block.zeroIdx},'SubSystem')
            % Subsystem that is not a BLOM block
            block.subsystem(block.zeroIdx) = true;
            block.reroute(block.zeroIdx) = true;
        elseif strcmp(block.blockType{block.zeroIdx},'From')
            % this is a from block
            block.fromBlock(block.zeroIdx) = true;       
            block.reroute(block.zeroIdx) = true;
        elseif strcmp(block.blockType{block.zeroIdx},'Demux')
            % DEMUX
            block.demux(block.zeroIdx) = true;
        elseif strcmp(block.blockType{block.zeroIdx},'Mux')
            % MUX
            block.mux(block.zeroIdx) = true;
        elseif strcmp(block.refBlock{block.zeroIdx}, 'BLOM_Lib/InputFromSimulink') ||...
                strcmp(block.refBlock{block.zeroIdx}, 'BLOM_Lib/InputFromWorkspace')
            % this is one of the BLOM input blocks
            block.inputBlock(block.zeroIdx) = true;
        elseif strcmp(block.refBlock{block.zeroIdx}, 'BLOM_Lib/ExternalFromSimulink') ||...
                strcmp(block.refBlock{block.zeroIdx}, 'BLOM_Lib/ExternalFromWorkspace')
            % this is one of the BLOM external blocks
            block.externalBlock(block.zeroIdx) = true;
        elseif strcmp(block.blockType{block.zeroIdx},'UnitDelay')
            % Unit Delay Label
            block.delay(block.zeroIdx) = true;
        elseif strcmp(block.blockType{block.zeroIdx},'Inport')
            parentOfBlock = get_param(block.handles(block.zeroIdx),'Parent');
            parentType = get_param(parentOfBlock,'BlockType');
            if strcmp(parentType,'SubSystem')
                block.subsystemInput(block.zeroIdx) = true;
                block.reroute(block.zeroIdx) = true;
            end
            
        else
            % store P and K matricies for the other blocks, as well as
            % boolean 1 if special function required and 0 ottherwise
            [P,K,specialFunPresence] = BLOM_Convert2Polyblock(currentBlockHandle);
            block.P{block.zeroIdx} = P;
            block.K{block.zeroIdx}= K;
            block.specialFunPresence(block.zeroIdx)=specialFunPresence;
        end

        if isempty(block.stepOutputIdx{block.zeroIdx})
            block.stepOutputIdx{block.zeroIdx} = zeros(length(sourcePorts.Outport),1);
        end
        

        % increase the index of block by one and populate currentBlockIndex
        currentBlockIndex = block.zeroIdx;
        block.zeroIdx = block.zeroIdx+1;
    else %case when the current outport's block is found but has not been added to block data
        currentBlockName = get_param(currentOutport,'Parent');
        currentBlockIndex = find(block.handles==currentBlockHandle);
        % FIX, find stuff for this case
    end
    
        
end

%% 
%======================================================================
%> @brief given list of outports and block, fill in the input field of
%> block
%>
%> More detailed description of the problem.
%>
%> @param block block structure
%> @param outportHandles all the outports found in BFS
%>
%> @retval block updated block structure with inputs field filled in
%======================================================================

function [block] = updateInputsField(block,stepVars,outportHandles)
    % get information for the block and port that this outport goes to
    for index = 1:length(outportHandles)
        currentOutport = outportHandles(index);
        % find the first index that corresponds to this outport in stepVars
        stepVarsIdx = find(stepVars.outportHandle==currentOutport,1);
             
        outportLine = get_param(currentOutport,'Line');
        destBlocks = get_param(outportLine,'DstBlockHandle');
        destPorts = get_param(outportLine,'DstPortHandle');
        % in the case of bounds and costs, we will not find block information
        % for it. Otherwise, this should always be true
        for i = 1:length(destBlocks)
            inportBlockIndex = block.handles==destBlocks(i);
            if any(inportBlockIndex)
                inportBlockPorts = get_param(block.handles(inportBlockIndex),'PortHandles');
                inportPorts = inportBlockPorts.Inport;
                inportIndex = inportPorts==destPorts(i);
                if isempty(block.stepInputIdx{inportBlockIndex})
                    block.stepInputIdx{inportBlockIndex} = zeros(length(inportPorts),1);
                end
                block.stepInputIdx{inportBlockIndex}(inportIndex) = stepVarsIdx;
            end
        end
    end

end


%%
%======================================================================
%> @brief given block and stepVars, create a large P and K matrix with the
%> indicies of stepVars.optVarIdx as the columns. 
%>
%> More detailed description of the problem.
%>
%> @param block block structure
%> @param stepVars stepVars structure
%>
%> @retval stepP P matrix of all variables. Columns of stepP correspond to
%> the optVarIdx field of stepVars
%> @retval stepK K matrix of all variables
%======================================================================

function [stepP,stepK] = combinePK(block,stepVars)
    stepP = [];
    stepK = [];
    optVarLength = max(stepVars.optVarIdx);
    
    for idx = 1:length(block.handles)
        % go through all the blocks and fill in what the P and K matrix
        % should be
        % idx = block index

        if ~isempty(block.P{idx})
            % only start filling in the P and K matrices if there are P and
            % K matrices to fill. not all blocks will have this field
            
            % find the current P and K matrix
            currentP = block.P{idx};
            currentK = block.K{idx};
            
            % get the full list of input and output indices
            inputsAndOutputsIdxs = [block.stepInputIdx{idx}; block.stepOutputIdx{idx}];

            % get the values of P and the corresponding columns and rows
            [rows,col,val] = find(block.P{idx});
            numRows = size(block.P{idx},1);
            % map the correct columns of P. First find which stepVars
            % indices the columns correspond to then find the optVarIdx
            % using the stepVars indices
            colIdx = inputsAndOutputsIdxs(col);
            optColIdx = stepVars.optVarIdx(colIdx);
            
            % make the newP with the proper column length and proper value
            % placement
            newP = sparse(rows,optColIdx,val,numRows,optVarLength);

            % concatenate the newP to the stepP
            stepP = vertcat(stepP,newP);
            stepK = blkdiag(stepK,block.K{idx});
            
        end
    end

    
end

%%
%================================================================
%> @brief updates block.inportIdxs and block.outportIdxs such that for
%> input and outport arrays, has a pointer to each scalar
%>
%> @param stepP the P matrix with all of columns as optVarIdx columns
%> @param stepK corresponding K matrix for stepP 
%>
%> @retval fullP the full P matrix for the entire problem for all time
%> steps
%> @retval fullK the full K matrix for the entire problem for all time
%> steps
%========================================================================

function [allP,allK] = createAllPK(stepP,stepK,stepVars,horizon,allVars)
    % get a mapping of which optVarIdx are relevant at each time step
    
    sizeTimes = length(stepVars.initTime);
    optInitTime = stepVars.initTime'*sparse(1:sizeTimes,stepVars.optVarIdx,1) > 0;
    optInterTime = stepVars.interTime'*sparse(1:sizeTimes,stepVars.optVarIdx,1) > 0;
    optFinalTime = stepVars.finalTime'*sparse(1:sizeTimes,stepVars.optVarIdx,1) > 0;
    
    % first get truncated P and K matrices with the relevant times and
    % variables
    [initP,initK] = trimPK(stepP,stepK,optInitTime);
    [interP,interK] = trimPK(stepP,stepK,optInterTime);
    [finalP,finalK] = trimPK(stepP,stepK,optFinalTime);
    
    % create the P and K matrices for all intermediate time steps
    if horizon>1
        interP_full = kron(speye(horizon-2),interP);
        interK_full = kron(speye(horizon-2),interK);

        fullP = blkdiag(initP,interP_full,finalP);
        allK = blkdiag(initK,interK_full,finalK);
    else
        fullP = initP;
        allK = initK;
    end
    
    % use allVars.PKOptVarIdxReroute to reroute columns of states
    [~,colPLength] = size(fullP);
    [~,~,rerouteIdx] = unique(allVars.PKOptVarIdxReroute);
    
    % create a matrix that finds which columns of the entire P (with
    % states) to add to the columns of P without states
    toAddCol = sparse(1:colPLength,rerouteIdx,ones(size(rerouteIdx)));
    
    % create the final allP that has proper states.
    allP = fullP*toAddCol;
end

%%
%================================================================
%> @brief updates block.inportIdxs and block.outportIdxs such that for
%> input and outport arrays, has a pointer to each scalar
%>
%> @param stepP the P matrix with all of columns as optVarIdx columns
%> @param stepK corresponding K matrix for stepP
%>
%> @retval trimP the P matrix for that time
%> @retval trimK K matrix for that time
%========================================================================

function [trimP,trimK] = trimPK(stepP,stepK,relevantTimes)
    % first create the P and K for initialTime
    trimP = stepP(:,relevantTimes);
    % remove initP rows if the row had a value in the columns that were
    % removed. 
    removeP_rows = any(stepP(:,~relevantTimes)~=0,2);
    trimP(removeP_rows,:) = [];
    % remove the corresponding columns of K
    trimK = stepK(:,~removeP_rows);
    % remove rows of K that had a value in the removed columns of K in the
    % previous line
    removeK_rows = any(stepK(:,removeP_rows)~=0,2);
    trimK(removeK_rows,:) = [];
end

%%
%================================================================
%> @brief updates block.inportIdxs and block.outportIdxs such that for
%> input and outport arrays, has a pointer to each scalar
%>
%> @param block 
%>
%> @retval block
%========================================================================

function block = expandInportOutportIdx(block)
    for blockIdx = 1:block.zeroIdx-1
        portDim = get_param(block.handles(blockIdx),'CompiledPortDimensions');
        inputsLength = sum(portDim.Inport(1:2:end).*portDim.Inport(2:2:end));
        outputsLength = sum(portDim.Outport(1:2:end).*portDim.Outport(2:2:end));
        
        if ~isempty(block.stepInputIdx{blockIdx})
            zeroIdx = 1;
            fullInputs = zeros(inputsLength,1);        
            for j = 1:length(block.stepInputIdx{blockIdx})
                currentInputLength = portDim.Inport(2*j-1).*portDim.Inport(2*j);
                fullInputs(zeroIdx:(zeroIdx+currentInputLength-1)) = ...
                    block.stepInputIdx{blockIdx}(j):(block.stepInputIdx{blockIdx}(j)+currentInputLength-1);
                zeroIdx = zeroIdx + currentInputLength;
            end
            block.stepInputIdx{blockIdx} = fullInputs;
        end
        
        
        if ~isempty(block.stepOutputIdx{blockIdx})
            zeroIdx = 1;
            fullOutputs = zeros(outputsLength,1);
            for j = 1:length(block.stepOutputIdx{blockIdx})
                currentOutputLength = portDim.Outport(2*j-1).*portDim.Outport(2*j);
                fullOutputs(zeroIdx:(zeroIdx+currentOutputLength-1)) = ...
                    block.stepOutputIdx{blockIdx}(j):(block.stepOutputIdx{blockIdx}(j)+currentOutputLength-1);
                zeroIdx = zeroIdx + currentOutputLength;
            end
            block.stepOutputIdx{blockIdx} = fullOutputs;
        end
        
        
    end
end


%%
%==========================================================================
%> @brief Reroutes pointers such that if System A points to System B, and 
%> System B points to System C, System A points to System C instead.
%> optVarIdx then reduced to remove gaps. 
%> e.g. [1;2;5;3;5]->[1;2;5;5;5]->[1;2;3;3;3].  
%> 
%> @param optVarIdx optVarIdx vector to be rerouted and reduced
%> 
%> @retval optVarIdx with optVarIdx rerouted and reduced
%==========================================================================
function optVarIdx = cleanupOptVarIdx(optVarIdx)

    for i = 1:length(optVarIdx)
        target = i;
        traversedTargets = zeros(size(optVarIdx));
        k = 1;
        while (target ~= optVarIdx(target))
            if(any(traversedTargets == target))

                traversedTargets(traversedTargets==0) = [];
                
                optVarIdx(traversedTargets) = min(traversedTargets);
                %error('ERROR: optVarIdx contains cycle')
                %slight redundancies, could be made faster with initial
                %check
            end
            traversedTargets(k) = target;
            k = k + 1;
            target = optVarIdx(target);
        end
        optVarIdx(i) = target;
    end
    [~,~,optVarIdx] = unique(optVarIdx);
end

%%
%====================================================================
%> @brief sets bounds in stepVars
%>
%> @param block all blocks from BFS
%> @param stepVars all stepVars from BFS
%>
%> @retval stepVars with fields initLowerBound, initUpperBound,
%> interLowerBound, interUpperBound, finalLowerBound, finalUpperBound
%> filled in
%======================================================================
function stepVars = setBounds(block, stepVars)
    
    searchIndices = find(block.bound);
    for currentBlockIndex = searchIndices'
        boundHandle = block.handles(currentBlockIndex);
        lowerBound = eval(get_param(boundHandle,'lb'));
        upperBound = eval(get_param(boundHandle,'ub'));
        
        % find out for which time steps these bounds are relevant
        init = strcmp(get_param(boundHandle, 'initial_step'), 'on');
        final = strcmp(get_param(boundHandle, 'final_step'), 'on');
        inter = strcmp(get_param(boundHandle, 'intermediate_step'), 'on');
        
        for varIdx = block.stepInputIdx{currentBlockIndex}'
            currentOutport = stepVars.outportHandle(varIdx);
            dimension = get_param(currentOutport,'CompiledPortDimensions');
            lengthOut = dimension(1)*dimension(2);
            % take into consideration current upper and lower bounds
            currentInitLB = stepVars.initLowerBound(varIdx);
            currentInterLB = stepVars.interLowerBound(varIdx);
            currentFinalLB = stepVars.finalLowerBound(varIdx);

            currentInitUB = stepVars.initUpperBound(varIdx);
            currentInterUB = stepVars.interUpperBound(varIdx);
            currentFinalUB = stepVars.finalUpperBound(varIdx);

            if init
                if lowerBound > currentInitLB
                    stepVars.initLowerBound(varIdx:(varIdx+lengthOut-1)) = lowerBound;
                end
                if upperBound < currentInitUB
                    stepVars.initUpperBound(varIdx:(varIdx+lengthOut-1)) = upperBound;
                end
            end

            if inter
                if lowerBound > currentInterLB
                    stepVars.interLowerBound(varIdx:(varIdx+lengthOut-1)) = lowerBound;
                end
                if upperBound < currentInterUB
                    stepVars.interUpperBound(varIdx:(varIdx+lengthOut-1)) = upperBound;
                end
            end

            if final
                if lowerBound > currentFinalLB
                    stepVars.finalLowerBound(varIdx:(varIdx+lengthOut-1)) = lowerBound;
                end
                if upperBound < currentFinalUB
                    stepVars.finalUpperBound(varIdx:(varIdx+lengthOut-1)) = upperBound;
                end
            end
        end
    end

end



%%
%====================================================================
%> @brief sets cost boolean in stepVars
%>
%> @param block all blocks from BFS
%> @param stepVars all stepVars from BFS
%>
%> @retval stepVars with fields initCost, interCost, finalCost filled in
%======================================================================
function stepVars = setCosts(block, stepVars)
    searchIndices = find(block.cost);
    for currentBlockIndex = searchIndices'
        costHandle = block.handles(currentBlockIndex);
        initCost =  strcmp(get_param(costHandle, 'initial_step'), 'on');
        interCost = strcmp(get_param(costHandle, 'intermediate_step'), 'on');
        finalCost = strcmp(get_param(costHandle, 'final_step'), 'on');   
        for varIdx = block.stepInputIdx{currentBlockIndex}'                   
            currentOutport = stepVars.outportHandle(varIdx);
            dimension = get_param(currentOutport,'CompiledPortDimensions');
            lengthOut = dimension(1)*dimension(2);
           
            stepVars.initCost(varIdx:(varIdx+lengthOut-1)) = initCost;
            stepVars.interCost(varIdx:(varIdx+lengthOut-1)) = interCost;
            stepVars.finalCost(varIdx:(varIdx+lengthOut-1)) = finalCost;
        end
    end
end
%%
%======================================================================
%> @brief traverse graph from cost and bound blocks then add fields to
%> stepVars to label relevance of variables at initial, intermediate, and 
%> final times.  For unit delay, all variables before the input are
%> relevant at all time steps before the relevance of the output variable
%>
%> @param block block structure
%>
%> @param stepVars all variables
%>
%> @retval stepVars stepVars with added boolean arrays as fields: initTime, 
%> interTime, finalTime
%=======================================================================

function stepVars = labelTimeRelevance(stepVars, block, inputAndExternalHandles)
    % the following stepVars fields state whether or not a variable is
    % relevant at that time
    stepVars.initTime = false(stepVars.zeroIdx-1,1);
    stepVars.interTime = false(stepVars.zeroIdx-1,1);
    stepVars.finalTime = false(stepVars.zeroIdx-1,1);
    
    blockIdxs = (1:block.zeroIdx-1)';
    
    startBlock = [blockIdxs(block.bound | block.cost); zeros(10,1)];
    startBlockIdx = 1;
    startBlockZeroIdx = sum(startBlock ~= 0)+1;
    while startBlockIdx ~= startBlockZeroIdx
        startBlockHandle = block.handles(startBlock(startBlockIdx));
        startBlockType = get_param(startBlockHandle, 'BlockType');

        if strcmp(startBlockType, 'UnitDelay')
            outputVarIdxs = block.stepOutputIdx{startBlock(startBlockIdx)};
            
            final = stepVars.finalTime(outputVarIdxs(1));
            inter = final || stepVars.interTime(outputVarIdxs(1));
            init = inter || stepVars.initTime(outputVarIdxs(1));
            
        else
            init = strcmp(get_param(startBlockHandle, 'initial_step'), 'on');
            final = strcmp(get_param(startBlockHandle, 'final_step'), 'on');
            inter = strcmp(get_param(startBlockHandle, 'intermediate_step'), 'on');
        end
        blocks = [startBlock(startBlockIdx); zeros(10,1)];
        
        zeroIdx = 2;
        idx = 1;
        while(blocks(idx)~=0)
            
            % if you've already been to current block previously or current
            % block is an input block (may be redundant), proceed to next loop iteration
            if(~any(blocks(idx)==blocks(1:idx-1)) && ~any(block.handles(blocks(idx))==inputAndExternalHandles))
                blockHandle = block.handles(blocks(idx));
                    
                if block.subsystem(blocks(idx))
                    block.names(idx) %should never be here
                elseif block.subsystemInput(blocks(idx))
                    % for the case of inport, we need to find the output
                    % connected to the subsystem and that specific inport
                    parentSubsystem = get_param(blockHandle,'Parent');
                    parentSubsystemPorts = get_param(parentSubsystem,'PortHandles');                    
                    portNum = eval(get_param(blockHandle,'Port'));
                    currentSubsystemInport = parentSubsystemPorts.Inport(portNum);
                    inportLine = get_param(currentSubsystemInport,'Line');
                    outportSource = get_param(inportLine,'SrcPortHandle');
                    inputs = find(stepVars.outportHandle == outportSource);

                elseif block.fromBlock(blocks(idx))
                    gotoTag = get_param(blockHandle, 'GotoTag');
                    name = get_param(blockHandle, 'Parent');
                    gotoBlock = find_system(name,'BlockType','Goto','GotoTag',gotoTag);
                    gotoPorts = get_param(gotoBlock{1},'PortHandles');
                    inputHandle = gotoPorts.Inport;
                    line = get_param(inputHandle, 'Line');
                    srcPortHandle = get_param(line, 'SrcPortHandle');
                    inputs = find(stepVars.outportHandle == srcPortHandle);
                    
                elseif block.delay(blocks(idx)) && (idx ~= 1)
                    inputVarIdxs = block.stepInputIdx{blocks(idx)};
                    inputFinal = stepVars.finalTime(inputVarIdxs(1));
                    inputInter = stepVars.interTime(inputVarIdxs(1));
                    inputInit = stepVars.initTime(inputVarIdxs(1));

                    if sum(startBlock(startBlockIdx:startBlockZeroIdx) == blocks(idx)) == 0 && ...
                            (inputFinal ~= final || inputInter ~= (inter || final) || inputInit ~= (final || inter || init))
                        startBlock(startBlockZeroIdx) = blocks(idx);
                        startBlockZeroIdx = startBlockZeroIdx + 1;
                        
                        if startBlockZeroIdx > length(startBlock)
                            startBlock = [startBlock zeros(1, length(startBlock))];
                        end
                    end
                    inputs = [];
                    
                else
                    inputs = block.stepInputIdx{blocks(idx)};
                end
                
                if(init)
                    stepVars.initTime(inputs) = true;
                end
                if(final)
                    stepVars.finalTime(inputs) = true;
                end
                if(inter)
                    stepVars.interTime(inputs) = true;
                end

                inputBlocks = stepVars.block(inputs);
                while any(block.subsystem(inputBlocks)) || any(block.mux(inputBlocks)) || any(block.demux(inputBlocks))
                    for inputIdx = 1:length(inputs)
                        if block.subsystem(inputBlocks(inputIdx))
                            subsystemHandle = block.handles(inputBlocks(inputIdx));
                            outportBlocks = find_system(subsystemHandle,'SearchDepth',1,'regexp','on','BlockType','Outport');
                            outportNum = stepVars.outportNum(inputs(inputIdx));
                            outportBlockPort = get_param(outportBlocks(outportNum), 'PortHandles');
                            outportInputHandle = outportBlockPort.Inport;
                            line = get_param(outportInputHandle, 'Line');
                            srcPortHandle = get_param(line, 'SrcPortHandle');
                            inputVector = find(stepVars.outportHandle == srcPortHandle);
                            vectorIndex = stepVars.outportIndex(inputs(inputIdx));
                            inputs(inputIdx) = inputVector(vectorIndex);
                        end
                        if block.mux(inputBlocks(inputIdx))
                            vectorIndex = find(block.stepOutputIdx{inputBlocks(inputIdx)}==inputs(inputIdx),1,'first');
                            inputVector = block.stepInputIdx{inputBlocks(inputIdx)};
                            inputs(inputIdx) = inputVector(vectorIndex);
                        end
                        if block.demux(inputBlocks(inputIdx))
                            vectorIndex = find(block.stepOutputIdx{inputBlocks(inputIdx)}==inputs(inputIdx),1,'first');
                            inputVector = block.stepInputIdx{inputBlocks(inputIdx)};
                            inputs(inputIdx) = inputVector(vectorIndex);
                        end
                        inputBlocks(inputIdx) = stepVars.block(inputs(inputIdx));  
                    end

                    if(init)
                        stepVars.initTime(inputs) = true;
                    end
                    if(final)
                        stepVars.finalTime(inputs) = true;
                    end
                    if(inter)
                        stepVars.interTime(inputs) = true;
                    end
                end
                
                %expand blocks when necessary
                if(length(blocks) < zeroIdx+length(inputs))
                    if(length(inputs) > length(blocks))
                        blocks = [blocks; zeros(2*length(inputs),1)];
                    else
                        blocks = [blocks; zeros(size(blocks))];
                    end
                end
               
                blocks(zeroIdx:zeroIdx+length(inputs)-1) = stepVars.block(inputs);
                zeroIdx = zeroIdx+length(inputs);
            end
            idx = idx+1;
        end
        
        startBlockIdx = startBlockIdx + 1;
    end

end



%%
%=====================================================================
%> @brief creates expands fields of block to account for multiple time
%> steps
%>
%> @param block completed block struct filled in for one time step
%>
%> @param timesteps number of timesteps
%>
%> @retval block struct with variables for all time steps
%=====================================================================
function block = expandBlock(block, horizon, stepVars, allVars)
    
    % adds a field to stepVars with row vectors containing pointers to
    % Allvars in column of appropriate timesteps
    % CHANGE ONLY APPLIES TO TEMP VARIABLE IN THIS FUNCTION CALL
    stepVars.allVarsIdxs = cell(stepVars.zeroIdx-1,1);
    [stepVars.allVarsIdxs{:}] = deal(zeros(horizon,1));
    for allVarsIdx = 1:allVars.totalLength
        stepVars.allVarsIdxs{allVars.stepVarIdx(allVarsIdx)}(allVars.timeStep(allVarsIdx)) = allVarsIdx;
    end
    
    block.allInputMatrix = cell(block.zeroIdx-1,1);
    block.allOutputMatrix = cell(block.zeroIdx-1,1);

    for blockIdx = 1:block.zeroIdx - 1
        block.allInputMatrix{blockIdx} = zeros(length(block.stepInputIdx{blockIdx}), horizon); 
        block.allOutputMatrix{blockIdx} = zeros(length(block.stepOutputIdx{blockIdx}), horizon); 
                
        block.allInputMatrix{blockIdx} = ...
            reshape([stepVars.allVarsIdxs{block.stepInputIdx{blockIdx}}], horizon, length(block.stepInputIdx{blockIdx}))';
        
        block.allOutputMatrix{blockIdx} = ...
            reshape([stepVars.allVarsIdxs{block.stepOutputIdx{blockIdx}}], horizon, length(block.stepOutputIdx{blockIdx}))';        
       
        block.allInputMatrix{blockIdx}(~any(block.allInputMatrix{blockIdx}~=0,2),:) = [];
        block.allOutputMatrix{blockIdx}(~any(block.allOutputMatrix{blockIdx}~=0,2),:) = [];
    end

end


%%
%=====================================================================
%> @brief creates allVars struct using separete time steps of stepVars
%>
%> @param stepVars completed stepVars filled in with time steps and bounds
%>
%> @retval allVars struct with variables separated by time steps
%=====================================================================
function allVars = createAllVars(stepVars,horizon)
%NOTE: currently fields of allVars are still tentative. Change to suit
%needs when necessary

    % first find length of allVars (this includes redundancies)
    initialLength = sum(stepVars.initTime);
    interLength = sum(stepVars.interTime);
    finalLength = sum(stepVars.finalTime);
    totalLength = initialLength+interLength*(horizon-2)+finalLength;
    
    allVars.totalLength = totalLength;
    
    % initialize allVars fields.
    allVars.lowerBound = zeros(totalLength,1);
    allVars.upperBound = zeros(totalLength,1);
    allVars.stepVarIdx = zeros(totalLength,1);
    allVars.optVarIdx = zeros(totalLength,1);
    allVars.timeStep = zeros(totalLength,1);
    allVars.cost = zeros(totalLength,1);
    
    % label what timestep each variable exists in
    allVars.timeStep(1:initialLength) = ones(initialLength,1);
    interSteps = ones(interLength,1) * (2:horizon-1);
    allVars.timeStep(initialLength+1:end-finalLength) = interSteps(:);
    allVars.timeStep(end-finalLength+1:end) = horizon * ones(finalLength,1);
    
    % set variable cost
    allVars.cost(1:initialLength) = stepVars.initCost(stepVars.initTime);
    allVars.cost(initialLength+1:end-finalLength) =...
        kron(ones(horizon-2,1),stepVars.interCost(stepVars.interTime));
    allVars.cost(end-finalLength+1:end) = stepVars.finalCost(stepVars.finalTime);   
    
    % because each variable has it's own time step, we can simply set
    % lowerBound and upperBound
    allVars.lowerBound(1:initialLength) = stepVars.initLowerBound(stepVars.initTime);
    allVars.lowerBound(initialLength+1:end-finalLength) =...
        kron(ones(horizon-2,1),stepVars.interLowerBound(stepVars.interTime));
    allVars.lowerBound(end-finalLength+1:end) = stepVars.finalLowerBound(stepVars.finalTime);
    
    allVars.upperBound(1:initialLength) = stepVars.initUpperBound(stepVars.initTime);
    allVars.upperBound(initialLength+1:end-finalLength) =...
        kron(ones(horizon-2,1),stepVars.interUpperBound(stepVars.interTime));
    allVars.upperBound(end-finalLength+1:end) = stepVars.finalUpperBound(stepVars.finalTime);
    
    % allVars.stepVarIdx points to the stepVarIdx that each index
    % corresponds to
    stepVarsIndices = 1:(stepVars.zeroIdx-1);
    allVars.stepVarIdx(1:initialLength) = stepVarsIndices(stepVars.initTime);
    allVars.stepVarIdx(initialLength+1:end-finalLength) =...
        kron(ones(horizon-2,1),stepVarsIndices(stepVars.interTime))';
    allVars.stepVarIdx(end-finalLength+1:end) = stepVarsIndices(stepVars.finalTime);

end

%%
%========================================================================
%> @brief creates allVars.optVarIdx and allVars.PKOptVarIdxReroute
%>
%> @param allVars to fill in
%>
%> @param block with block.stepOutputIdx and allInputMatrix already filled
%> in.  currently done in expandBlock
%>
%> @param stepVars with initTime, interTime, finaltime filled in
%> 
%> @param horizon for problem
%>
%> @retval allVars with allVars.optVarIdx and allVars.PKOptVarIdxReroute
%> filled in
%======================================================================
function allVars = allOptVarIdxs(allVars,block,stepVars,horizon)
    
    initialLength = sum(stepVars.initTime);
    interLength = sum(stepVars.interTime);
    finalLength = sum(stepVars.finalTime);
    totalLength = initialLength+interLength*(horizon-2)+finalLength;

    % init Time
    stepVarLength = stepVars.zeroIdx - 1;
    allVars.optVarIdx(1:initialLength) = stepVars.optVarIdx(stepVars.initTime);
    % inter Time
    repeatedInterOptVarIdx = repmat(stepVars.optVarIdx(stepVars.interTime),horizon-2,1);
    incrementInterOptVarIdxMat = stepVarLength * ones(interLength,1) * (1:horizon-2);
    incrementInterOptVarIdxVector = incrementInterOptVarIdxMat(:);
    allVars.optVarIdx(initialLength+1:initialLength+interLength*(horizon-2)) = ...
        repeatedInterOptVarIdx + incrementInterOptVarIdxVector;
    % final Time
    allVars.optVarIdx(end-finalLength+1:end) = stepVars.optVarIdx(stepVars.finalTime) + stepVarLength*(horizon-1);
        
    % reroute optVarIdx for delay blocks such that output of 1/z at timestep
    % k+1 is same as optvaridx of inut at timestep k
    [~,~,allVars.optVarIdx] = unique(allVars.optVarIdx);
        
    allVars.PKOptVarIdxReroute = (1:max(allVars.optVarIdx))';
    for idx = initialLength+1:totalLength
        stepVarIdx = allVars.stepVarIdx(idx);
        if stepVars.state(stepVarIdx)
            oldOptVarIdx = allVars.optVarIdx(idx);
            blockIdx = stepVars.block(stepVarIdx);
            timeStep = allVars.timeStep(idx);
            blockInputOutputIdx = find(block.stepOutputIdx{blockIdx}== stepVarIdx,1);
            newOptVarIdx = allVars.optVarIdx(block.allInputMatrix{blockIdx}(blockInputOutputIdx,timeStep-1));

            allVars.optVarIdx(allVars.optVarIdx == oldOptVarIdx) = newOptVarIdx;
            allVars.PKOptVarIdxReroute(oldOptVarIdx) = newOptVarIdx;
        end
    end
    
    [~,~,allVars.optVarIdx] = unique(allVars.optVarIdx);    

end


%%
%======================================================================
%> @brief Convert BLOM 2.0 variables into a structure readable by BLOM 1.0
%>
%> @param block: block structure
%> @param stepVars: all variables
%>
%> @retval ModelSpec structure that will feed into BLOM 2.0
%=======================================================================


% NOTE, REPLACE VARARGIN WITH ACTUAL VARIABLES. THIS IS JUST A PLACEHOLDER
% FOR NOW
function [ModelSpec] = convert2ModelSpec(name,horizon,integ_method,dt,options,stepVars,allVars,block,allP,allK,varargin)
    ModelSpec.name = name;
    ModelSpec.horizon = horizon;
    ModelSpec.integ_method = integ_method;
    ModelSpec.dt = dt;
    ModelSpec.options = options;

    numOptVars = max(allVars.optVarIdx);

    
    ModelSpec.in_vars = false(numOptVars,1);
%     ModelSpec.in_vars(allVars.optVarIdx) = stepVars.input(allVars.stepVarIdx);
    ModelSpec.ex_vars = false(numOptVars,1);
%     ModelSpec.ex_vars(allVars.optVarIdx) = stepVars.external(allVars.stepVarIdx);
    for ii = 1:allVars.totalLength
       ModelSpec.in_vars(allVars.optVarIdx(ii)) = ModelSpec.in_vars(allVars.optVarIdx(ii)) || stepVars.input(allVars.stepVarIdx(ii));
       ModelSpec.ex_vars(allVars.optVarIdx(ii)) = ModelSpec.ex_vars(allVars.optVarIdx(ii)) || stepVars.external(allVars.stepVarIdx(ii));
    end

    
    ModelSpec.AAs = {allP};
    ModelSpec.Cs = {allK};
    
    varNameParentChar = char(block.names(stepVars.block(allVars.stepVarIdx)));
    varNameParent = cellstr(varNameParentChar(:,(length(name)+2):end));
    for ii = 1:length(varNameParent)
        varNameParent{ii}(varNameParent{ii}=='/') = '_';
    end
    stepVarOutputNum = zeros(stepVars.zeroIdx-1,1);
    for blockIdx = 1:block.zeroIdx-1      
        for outputIdx = 1:length(block.stepOutputIdx{blockIdx})
            stepVarOutputNum(block.stepOutputIdx{blockIdx}(outputIdx)) = outputIdx;
        end
    end
    varNameOutputNum = num2str(stepVarOutputNum(allVars.stepVarIdx));
    varNameTimeStep = num2str(allVars.timeStep);
    varNamePortNum = num2str(stepVars.outportNum(allVars.stepVarIdx));
    varNameVecIdx = num2str(stepVars.outportIndex(allVars.stepVarIdx));
    varNameAllVars = strcat('BL_',varNameParent, '.Out', varNameOutputNum, '.t', varNameTimeStep,'.port', varNamePortNum  ,'.vecIdx', varNameVecIdx, ';');
    
    
    %create all_names field
    ModelSpec.all_names = cell(numOptVars,1);
    for idx = 1:allVars.totalLength
        ModelSpec.all_names{allVars.optVarIdx(idx)} = ...
            strcat(ModelSpec.all_names{allVars.optVarIdx(idx)}, varNameAllVars{idx});
    end
    
    for idx = 1:size(ModelSpec.all_names)
       ModelSpec.all_names{idx} =  ModelSpec.all_names{idx}(1:end-1);
    end
   
    


    %create all_names_struct
    num_terms = cellfun(@length, strfind(ModelSpec.all_names,';')) + 1; % number of ';'
    terms_so_far = [0; cumsum(num_terms)]; % is number of multiple names
    if isempty(ModelSpec.all_names)
       all_fields = cell(1,5);
    else
        all_fields = textscan([ModelSpec.all_names{:}],'BL_%sOut%dt%dport%dvecIdx%d','Delimiter','.;');
    end
    
    vec_idx = zeros(terms_so_far(end),1); % preallocate vec_idx
    vec_idx(terms_so_far(1:end-1)+1) = 1:numOptVars; % first of each
    twoterms = find(num_terms == 2);
    vec_idx(terms_so_far(twoterms)+2) = twoterms; % 2nd of each
    multiterms = find(num_terms > 2); % multiples, should be fewer of these
    for i = 1:length(multiterms)
        vec_idx(terms_so_far(multiterms(i))+2 : ...
            terms_so_far(multiterms(i)+1)) = multiterms(i);
    end
    
    ModelSpec.all_names_struct.terms_so_far = terms_so_far;
    ModelSpec.all_names_struct.all_fields = all_fields;
    ModelSpec.all_names_struct.vec_idx = vec_idx; 

    
    %create inequality polyblocks
    optVarsBounds = zeros(numOptVars,2);   
    optVarsBounds(allVars.optVarIdx,1) = allVars.lowerBound;
    optVarsBounds(allVars.optVarIdx,2) = allVars.upperBound;
    lowBounded = ~isinf(optVarsBounds(:,1));
    upBounded = ~isinf(optVarsBounds(:,2));
    numBounds = sum(lowBounded)+sum(upBounded);
    lowBounds = optVarsBounds(lowBounded,1);
    upBounds = optVarsBounds(upBounded,2); 
    
    if numBounds > 0
        ModelSpec.ineq.AAs = {sparse(1:numBounds, [find(upBounded)' find(lowBounded)'], ones(1,numBounds), numBounds+1, numOptVars)};
        
        ineqCsM = [1:numBounds 1:numBounds];
        ineqCsN = [1:numBounds (numBounds+1)*ones(1,numBounds)];
        ineqCsVals = [ones(1,sum(upBounded)) -ones(1,sum(lowBounded)) -upBounds' lowBounds'];
        ModelSpec.ineq.Cs = {sparse(ineqCsM, ineqCsN, ineqCsVals,numBounds, numBounds+1)};
    else
        % FIX: check if this works with BLOM1 solver interfaces for
        % problems without inequalities
        ModelSpec.ineq.AAs = {[]};
        ModelSpec.ineq.Cs = {[]};
    end
    
    %create cost polyblocks
    optVarCosts = zeros(numOptVars,1);
    optVarCosts(allVars.optVarIdx) = allVars.cost;
    nonzeroOptVarCosts = optVarCosts~=0;
    if any(nonzeroOptVarCosts)
        ModelSpec.cost.A = sparse(1:sum(nonzeroOptVarCosts), find(nonzeroOptVarCosts),...
            optVarCosts(nonzeroOptVarCosts), sum(nonzeroOptVarCosts), numOptVars);
        ModelSpec.cost.C = ones(1, sum(nonzeroOptVarCosts));
    elseif numOptVars > 0
        % cost is 0, this is a feasibility problem
        ModelSpec.cost.A = sparse(1, numOptVars);
        ModelSpec.cost.C = 0;
    else
        % no optimization variables?
        ModelSpec.cost.A = [];
        ModelSpec.cost.C = [];
    end
    
    %create A,C polyblocks
    ModelSpec.A = [ModelSpec.cost.A; ModelSpec.ineq.AAs{1}; ModelSpec.AAs{1}];
    ModelSpec.C = blkdiag(ModelSpec.cost.C, ModelSpec.ineq.Cs{1}, ModelSpec.Cs{1});
    [len_cost_A,~] = size(ModelSpec.cost.A);
    [len_ineq_A,~] = size(ModelSpec.ineq.AAs{1});
    [len_eq_A,~] = size(ModelSpec.AAs{1});
    [len_cost_C,~] = size(ModelSpec.cost.C);
    [len_ineq_C,~] = size(ModelSpec.ineq.Cs{1});
    [len_eq_C,~] = size(ModelSpec.Cs{1});        
    
    ModelSpec.ineq_start_A = len_cost_A + 1;
    ModelSpec.ineq_end_A = len_cost_A + len_ineq_A;
    ModelSpec.eq_start_A = len_cost_A + len_ineq_A + 1;
    ModelSpec.eq_end_A = len_cost_A + len_ineq_A + len_eq_A;
    ModelSpec.ineq_start_C = len_cost_C + 1;
    ModelSpec.ineq_end_C = len_cost_C + len_ineq_C;
    ModelSpec.eq_start_C = len_cost_C + len_ineq_C + 1;
    ModelSpec.eq_end_C = len_cost_C + len_ineq_C + len_eq_C;
    
    %create all_state_vars
    % FIX: this should probably not mix stepVars and allVars, since some
    % of stepVars might not be included in allVars for first time step
    ModelSpec.all_state_vars = sparse((allVars.optVarIdx(stepVars.state)), ones(sum(stepVars.state),1), ones(sum(stepVars.state),1), numOptVars,1);
 
end
