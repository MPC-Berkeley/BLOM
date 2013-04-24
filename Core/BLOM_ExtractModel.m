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

function [ModelSpec,block,stepVars,allVars] = BLOM_ExtractModel(name,horizon,dt,integ_method,options)
    % load system. does nothing if model is not open
    load_system(name);
    % evaluate model to get dimensions
    eval([name '([],[],[],''compile'');']); 
    try
        [boundHandles,costHandles,inputAndExternalHandles] = findBlocks(name);
        [block,stepVars,stop] = ...
            searchSources(boundHandles,costHandles,inputAndExternalHandles,name);
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
        % [allP,allK] = createAllPK(stepP,stepK,stepVars,horizon,allVars);

        % create initial fields for ModelSpec
        ModelSpec.Name = name;
        ModelSpec.integ_method = integ_method;
        ModelSpec.dt = dt;
        ModelSpec.horizon = horizon;
        ModelSpec.options = options;
        
        % convert to ModelSpec - this part will merge BLOM 2.0 and BLOM 1.0
        [ModelSpec] = convert2ModelSpec(ModelSpec,stepVars,block);
        
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
                inportsOutport(index) = get_param(currentLine,'SrcPortHandle');
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

        % just a placeholder for ModelSpec so that MATLAB does not complain
        ModelSpec = 1;
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
    stepVars.outportIndex = zeros(initialSize,1); % index of specific outport. normally just 1 (will be more than 1 if the outport is a vector)
    stepVars.optVarIdx = (1:initialSize)'; % points to the true optimization variable index. (will usually point to itself, otherwise, some "true" variable)
    stepVars.cost = zeros(initialSize,1); %default cost is zero. 1 if we see that it's part of the cost
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
    block.dimensions = cell(initialSize,1); % dimensions of each outport. first value is outport #, then second two values are dimensions of outport
    block.bound = false(initialSize,1); % indicator to whether block is a bound block
    block.cost = false(initialSize,1); % indicator to whether block is a cost block
    
    % find all lines connected to costs and bounds and then get outport
    % ports from there. fill in stepVars and block structures as necessary
    
    % iZero is index of first zero of outportHandles
    iZero = 1;
    outportHandles = zeros(initialSize,1);
    
    % get all outports connected to the bounds 
    state = 'bound';
    for i = 1:length(boundHandles)
        portH = get_param(boundHandles(i),'PortHandles');
        % costs and bounds should only have one inport and line
        currentInport = [portH.Inport];
        [block, ~] = updateBlock(block,currentInport);
        [outportHandles,iZero,stepVars,block] = ...
            updateVars(currentInport,outportHandles,iZero,...
            stepVars,block,state,iZero,portH,boundHandles(i));
    end
 
    
    % get all outports connected to costs
    state = 'cost';
    for i = 1:length(costHandles)
        portH = get_param(costHandles(i),'PortHandles');
        % costs and bounds should only have one inport and line
        currentInports = [portH.Inport];
        [block, ~] = updateBlock(block, currentInports);
        [outportHandles,iZero,stepVars,block] = ...
            updateVars(currentInports,outportHandles,iZero,...
            stepVars,block,state);
    end
    
    
    % iOut is current outport handle we are looking at
	iOut = 1;
    while 1
        if outportHandles(iOut) == 0
            % all handles have been found and searched
            break
        elseif any(outportHandles(iOut)==sOutportHandles)
            % if the current handle is equal to any of the external or input
            % from simulink, then continue with the loop and do not look
            % for more branches from that block
            state = 'fromSimulink';
            sourceInports = 0;
            sourcePorts = 0;
            [outportHandles,iZero,stepVars,block] = ...
                    updateVars(sourceInports,outportHandles,iZero,...
                    stepVars,block,state,iOut,sourcePorts);
            
            iOut = iOut+1;
            continue
        elseif iOut > length(outportHandles)
            % should not reach here, just in case
            break
        elseif ~outportHandles(end)==0
            % no more space for more handles, allocate more space
            outportHandles = [outportHandles; zeros(length(outportHandles)*2,1)];
        end
        
        sourceBlock = get_param(outportHandles(iOut),'Parent');
        sourcePorts = get_param(sourceBlock,'PortHandles');
        sourceInports = [sourcePorts.Inport];
        sourceType = get_param(sourceBlock,'BlockType');
        refBlock = get_param(sourceBlock,'ReferenceBlock');

        
        if strcmp(sourceType,'SubSystem')
            if isempty(refBlock)
                % if the current block is a subsystem and not from BLOM, 
                % want to look under the subsystem and get the appropriate
                % blocks there
                
                state = 'intoSubsys';
                sourceOutports = [sourcePorts.Outport];
                [outportHandles,iZero,stepVars,block] = ...
                    updateVars(sourceOutports,outportHandles,iZero,...
                    stepVars,block,state,iOut,sourcePorts);

            elseif strncmp(refBlock,'BLOM_Lib',8)
                % subsystem that's a BLOM block
                state = 'BLOM';
                [outportHandles,iZero,stepVars,block] = ...
                    updateVars(sourceInports,outportHandles,iZero,...
                    stepVars,block,state,iOut,sourcePorts);
            end
                
        elseif strcmp(sourceType,'From')
            % the current block is a from block, find goto block and see
            % what is connected to the goto block
            % NOTE: sourceInports is really the tag of the from block ONLY
            % in this From block case. in the updateVars code, it takes
            % this into consideration and finds the proper blocks
            sourceInports{1} = get_param(sourceBlock,'GotoTag');
            sourceInports{2} = name;

            state = 'from';
            [outportHandles,iZero,stepVars,block] = ...
                updateVars(sourceInports,outportHandles,iZero,...
                stepVars,block,state,iOut,sourcePorts);
            
        elseif strcmp(sourceType,'Demux')
            % block is a demux. Here we want to fill in field sameOptVar to
            % point to the original variable
            
            state = 'demux';
            sourceOutports = [sourcePorts.Outport];
            [outportHandles,iZero,stepVars,block] = ...
                updateVars(sourceInports,outportHandles,iZero,...
                stepVars,block,state,iOut,sourcePorts,sourceOutports);
            
        elseif strcmp(sourceType,'Mux')
            % block is a mux. Here we want to fill in field sameOptVar to
            % point to the original variable
            
            state='mux';
            [outportHandles,iZero,stepVars,block] = ...
                updateVars(sourceInports,outportHandles,iZero,...
                stepVars,block,state,iOut,sourcePorts);
            
        elseif strcmp(sourceType,'UnitDelay')       
            state = 'unitDelay';
            [outportHandles,iZero,stepVars,block] = ...
                updateVars(sourceInports,outportHandles,iZero,...
                stepVars,block,state,iOut,sourcePorts);
            
        elseif length(sourceInports) > 1
            % currently do nothing special if there is more than one input
            % except look for what's connected. 
            % FIX: check to see which inports affect which outports. Not
            % all inports may be relevant
            state = 'norm';
            [outportHandles,iZero,stepVars,block] = ...
                updateVars(sourceInports,outportHandles,iZero,...
                stepVars,block,state,iOut,sourcePorts);
            
        elseif isempty(sourceInports)
            % if there are no inports, no need to search this outport
            % anymore. However, if the block is an inport of a subsystem,
            % we must see what is connected to that block
            parentOfBlock = get_param(sourceBlock,'Parent');
            if ~strcmp(parentOfBlock,name)
                % FIX: Is this still necessary?
                parentType = get_param(parentOfBlock,'BlockType');
                if strcmp(sourceType,'Inport') && strcmp(parentType,'SubSystem')
                    % in this case, there may actually be more inports
                    % connected
                    sourcePorts = get_param(parentOfBlock,'PortHandles');
                    sourceInports = [sourcePorts.Inport];
                    % find the specific inport index of 
                    portNumber = eval(get_param(sourceBlock,'Port'));
                    % storing the inport itself is not important. what we
                    % need to store is what is connected to the inport of
                    % the subsystem. so, store information for that inport,
                    % but continue BFS from the subsystem itself.

                    if ~isempty(sourceInports)
                        state = 'subSysInport';
                        [outportHandles,iZero,stepVars,block] = ...
                            updateVars(sourceInports(portNumber),outportHandles,iZero,...
                            stepVars,block,state,iOut,sourcePorts);
                    else
                        % This case should NOT be reached because of an
                        % inport exists, there is almost certainly a
                        % corresponding inport for the entire subsystem
                        % itself. It is here as asafeguard. 
                        stepVarsState = 'normal';
                        fprintf('Have you forgotten to connect one of the inports to subsystem: %s? ', parentOfBlock);
                    end
                else
                    stepVarsState = 'normal';
                    [block,currentBlockIndex] =...
                        updateBlock(block,outportHandles(iOut));
                    [stepVars,block] = updateStepVars(stepVars,...
                        block,currentBlockIndex,outportHandles(iOut),stepVarsState);
                end
            else
                stepVarsState = 'normal';
                [block,currentBlockIndex] =...
                    updateBlock(block,outportHandles(iOut));
                [stepVars,block] = updateStepVars(stepVars,...
                    block,currentBlockIndex,outportHandles(iOut),stepVarsState);
                iOut = iOut+1;
                continue
            end
        else
            % in the case of one inport, it must affect the outport, so
            % find the relevant outports
            [outportHandles,iZero,stepVars,block] = ...
                updateVars(sourceInports,outportHandles,iZero,...
                stepVars,block,2,iOut,sourcePorts);
        end
        % FIX: should check to see if BLOM supports the blocks that is
        % found right here
        iOut = iOut+1;
    end
    
    
    if any(outportHandles==0)
        outportHandles = outportHandles(1:iZero-1);
    end
    
    block = updateInputsField(block,stepVars,outportHandles);
    
    %check if any handles are -1 and remove them
    if any(outportHandles==-1)
        % FIX PRINT STATEMENT. NOT NECESSARILY TRUE
        fprintf('One or more of the blocks are missing a connection\n');
        stop = 1;
    end
    
    % remove empty and 0 entries in blocks
    for field={'names', 'P','K','stepInputIdx','stepOutputIdx','dimensions'}
        block.(field{1}) = block.(field{1})(1:(block.zeroIdx-1));
    end
    block.handles = block.handles(1:(block.zeroIdx-1));
    block.bound = block.bound(1:(block.zeroIdx-1));
    block.cost = block.cost(1:(block.zeroIdx-1));

    
    % remove empty and 0 entries in stepVars
    for field={'block', 'outportNum','outportHandle','outportIndex',...
        'optVarIdx','cost','state','input','external','initLowerBound','initUpperBound',...
        'interLowerBound','interUpperBound','finalLowerBound','finalUpperBound'}
        stepVars.(field{1}) = stepVars.(field{1})(1:(stepVars.zeroIdx-1));
    end

    % FIX: need to find some way to remove all -1 handles. using setdiff
    % with [-1] reorders all the outport handles and puts it in ascending
    % order
    % outportHandles = setdiff(outportHandles,[-1]);
    
end

%======================================================================
%> @brief From given inports, see which outports are relevant. Also, given
%> the current outport, populate the block and stepVars structures
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

function [outportHandles,iZero,stepVars,block] =...
    updateVars(inports,existingOutports,iZero,stepVars,...
    block,state,varargin)

    outportHandles = existingOutports;

    
    if ~isempty(varargin)
        iOut = varargin{1};
        sourcePorts = varargin{2};
        currentOutport = existingOutports(iOut);
    end
    
    if strcmp(state,'intoSubsys') %FIX: may want to look into special BLOM blocks here
        % if there's a subsystem, inports is actually an array of the
        % outports of subsystem
        parent = get_param(currentOutport,'Parent');
        index = inports==currentOutport;
        outportBlocks = find_system(parent,'SearchDepth',1,'regexp','on','BlockType','Outport');
        handle = get_param(outportBlocks{index},'Handle');
        portH = get_param(handle,'PortHandles');
        inports = [portH.Inport];
    end
    
    if strcmp(state,'from')
        tag = inports{1};
        name = inports{2};
        % there should only be one goto block associated with this from
        % block
        gotoBlock = find_system(name,'BlockType','Goto','GotoTag',tag);
        gotoPorts = get_param(gotoBlock{1},'PortHandles');
        % note: there should only be one outport associated with each goto
        % block
        inports = [gotoPorts.Inport];
        
    end
   
    % found outports connected to inports provided
    if inports ~= 0
        allOutportsFound=zeros(length(inports),1);
        for i = 1:length(inports);
            currentLine = get_param(inports(i),'Line');
            % this gives the all the outports connected to this line
            outportFound = get_param(currentLine,'SrcPorthandle');
            allOutportsFound(i) = outportFound;
            outLength = length(outportFound);
            % in case outportHandles is too short 
            if outLength > length(outportHandles)-iZero+1;
                if outlength > length(outportHandles)
                    outportHandles = [outportHandles; zeros(outLength*2,1)];
                else
                    outportHandles = [outportHandles; ...
                        zeros(length(outportHandles),1)];
                end
            end

            diffOutports = setdiff(outportFound,outportHandles);
            diffLength = length(diffOutports);
            outportHandles(iZero:(diffLength+iZero-1)) = diffOutports;
            iZero = iZero + diffLength;
        end
    end
    
    switch state
        case 'bound'
            % bound case. for all the outports found here, include the
            % bounds for these
            stepVarsState = 'bound';
            boundHandle = varargin{3};
            
            % go through each of the outports connected through bound.
            % Because bound has only one inport, the for loop above only
            % goes through one iteration, so the variable outportFound
            % will include all the outports found
            for currentOutport = outportFound
                [block,currentBlockIndex] =...
                    updateBlock(block,currentOutport);
                [stepVars,block] = updateStepVars(stepVars,...
                    block,currentBlockIndex,currentOutport,stepVarsState,boundHandle);
            end
            
        case 'cost'
            stepVarsState = 'cost';
            % go through each of the outports connected through cost
            % Because cost has only one inport, the for loop above only
            % goes through one iteration, so the variable outportFound
            % will include all the outports found
            for currentOutport = outportFound
                [block,currentBlockIndex] =...
                    updateBlock(block,currentOutport);
                [stepVars,block] = updateStepVars(stepVars,...
                    block,currentBlockIndex,currentOutport,stepVarsState);
            end
            
        case 'fromSimulink'
            % for the fromSimulink blocks, we need a special case to look
            % at that specific outport even though we stop the BFS here
            
%             stepVarsState = 'normal';
%             currentOutport = existingOutports(iOut);
%             [block,currentBlockIndex] =...
%                 updateBlock(block,currentOutport);
%             [stepVars,block] = updateStepVars(stepVars,...
%                 block,currentBlockIndex,currentOutport,stepVarsState);
%             
            sourceBlock = get_param(outportHandles(iOut),'Parent');
            refBlock = get_param(sourceBlock,'ReferenceBlock');

            if strcmp(refBlock, 'BLOM_Lib/InputFromSimulink') || ...
               strcmp(refBlock, 'BLOM_Lib/InputFromWorkspace')
                stepVarsState = 'input';

                currentOutport = existingOutports(iOut);
                [block,currentBlockIndex] =...
                    updateBlock(block,currentOutport);
                [stepVars,block] = updateStepVars(stepVars,...
                    block,currentBlockIndex,currentOutport,stepVarsState);
            elseif strcmp(refBlock, 'BLOM_Lib/ExternalFromWorkspace') || ...
                   strcmp(refBlock, 'BLOM_Lib/ExternalFromSimulink')
            
                stepVarsState = 'external';

                currentOutport = existingOutports(iOut);
                [block,currentBlockIndex] =...
                    updateBlock(block,currentOutport);
                [stepVars,block] = updateStepVars(stepVars,...
                    block,currentBlockIndex,currentOutport,stepVarsState);
            else
                %should never reach here
            end
            
        case 'subSysInport'
            % Here we fill in the optVarIdx field have the inport point to
            % the outport connected to the inport of the subsystem

            % the current outport changes here because we want the outport
            % connected to the subsystem which is connected to the inport
            % to the subsystem to have a reference that it's a duplicate
            % variable
            stepVarsState = 'rememberIndex';
            currentOutport = outportFound;
            [block,currentBlockIndex] =...
                updateBlock(block,currentOutport);
            [stepVars,block,sameOptIndex] = updateStepVars(stepVars,...
                block,currentBlockIndex,currentOutport,stepVarsState);
            
            stepVarsState = 'subSysInport';
            currentOutport = existingOutports(iOut);
            [block,currentBlockIndex] =...
                updateBlock(block,currentOutport);
            [stepVars,block] = updateStepVars(stepVars,...
                    block,currentBlockIndex,currentOutport,stepVarsState,sameOptIndex);
        case 'from'
            % in this case, we want to be able to fill in the sameOptVar
            % field in stepVars so that we don't create a separate
            % optimization variable but rather know that this is a duplicate
            stepVarsState = 'rememberIndex';
            % the currentOutport in this case is the outport that goes to
            % the GoTo block
            currentOutport = outportFound;
            [block,currentBlockIndex] =...
                updateBlock(block,currentOutport);
            [stepVars,block,sameOptIndex] = updateStepVars(stepVars,...
                    block,currentBlockIndex,currentOutport,stepVarsState);
            % save information for from block here
            stepVarsState = 'from';
            currentOutport = existingOutports(iOut);
            [block,currentBlockIndex] =...
                updateBlock(block,currentOutport);
            [stepVars,block] = updateStepVars(stepVars,...
                    block,currentBlockIndex,currentOutport,stepVarsState,sameOptIndex);
                
        case 'intoSubsys'
            % for this case, we want to fill in the sameOptVar variable for
            % the subsystems outport and say the original variable is the
            % outport of whatever is inside the subsystem
            stepVarsState = 'rememberIndex';
            % the currentOutport in this case is the outport that goes to
            % the subsystem
            currentOutport = outportFound;
            [block,currentBlockIndex] =...
                updateBlock(block,currentOutport);
            [stepVars,block,sameOptIndex] = updateStepVars(stepVars,...
                block,currentBlockIndex,currentOutport,stepVarsState);
            
            % save information for subsystem outport here
            stepVarsState = 'intoSubsys';
            currentOutport = existingOutports(iOut);
            [block,currentBlockIndex] =...
                updateBlock(block,currentOutport);
            [stepVars,block] = updateStepVars(stepVars,...
                block,currentBlockIndex,currentOutport,stepVarsState,sameOptIndex);
                
        case 'demux'
            % the original variable is going to be the outport connected to
            % the demux
            
            stepVarsState = 'rememberIndex';
            % should only be one outport connected to the demux.
            currentOutport = outportFound;
            [block,currentBlockIndex] =...
                updateBlock(block,currentOutport);
            [stepVars,block,sameOptIndex] = updateStepVars(stepVars,...
                block,currentBlockIndex,currentOutport,stepVarsState);
            
            % fill in optVarIdx for the demux's outports to the outport
            % that goes into demux
            sourceOutports = varargin{3};
            for i = 1:length(sourceOutports)
                % store information for the rest of the demux outports
                stepVarsState = 'demux';
                [block,currentBlockIndex] =...
                    updateBlock(block,sourceOutports(i));
                [stepVars,block] = updateStepVars(stepVars,...
                    block,currentBlockIndex,sourceOutports(i),stepVarsState,sameOptIndex,i);
            end
                
            
            
        case 'mux'           
            % the original variables are the outports connected to the mux
            mux_optVarIdx = zeros(length(allOutportsFound),1);
            for idx=1:length(allOutportsFound)
                   currentOutport = allOutportsFound(idx);                
                   stepVarsState = 'rememberIndex';
                    [block,currentBlockIndex] =...
                        updateBlock(block,currentOutport);
                    [stepVars,block,sameOptIndex] = updateStepVars(stepVars,...
                        block,currentBlockIndex,currentOutport,stepVarsState);
                    mux_optVarIdx(idx) = sameOptIndex;
             end    
              % save information for mux outport here   
                stepVarsState = 'mux';
                currentOutport=existingOutports(iOut);
                [block,currentBlockIndex] =...
                    updateBlock(block,currentOutport);
                [stepVars,block] = updateStepVars(stepVars,...
                    block,currentBlockIndex,currentOutport,stepVarsState,mux_optVarIdx);  

        case 'unitDelay'
            stepVarsState = 'unitDelay';
            %update block structure
            [block,currentBlockIndex] =...
                updateBlock(block,currentOutport);
            [stepVars,block] = updateStepVars(stepVars,...
                block,currentBlockIndex,currentOutport,stepVarsState);
            
        otherwise 
            % not looking at bounds or costs
            stepVarsState = 'normal';
            %update block structure
            [block,currentBlockIndex] =...
                updateBlock(block,currentOutport);
            [stepVars,block] = updateStepVars(stepVars,...
                block,currentBlockIndex,currentOutport,stepVarsState);
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
    block,currentBlockIndex,currentOutport,stepVarsState,varargin)
%this function populates stepVars structure
%state can be 'bound', 'cost' 
    if sum(any(stepVars.outportHandle==currentOutport)) == 0
        %only add outports which haven't been looked at
        dimension = get_param(currentOutport,'CompiledPortDimensions');
        lengthOut = dimension(1)*dimension(2);
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
                        'cost'}
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
        
        switch stepVarsState
            case 'bound'
                % fill in bound parameters for different time steps
                boundHandle = varargin{1};
                lowerBound = eval(get_param(boundHandle,'lb'));
                upperBound = eval(get_param(boundHandle,'ub'));

                % find out for which time steps these bounds are relevant
                init = strcmp(get_param(boundHandle, 'initial_step'), 'on');
                final = strcmp(get_param(boundHandle, 'final_step'), 'on');
                inter = strcmp(get_param(boundHandle, 'intermediate_step'), 'on');
                
                if init
                    stepVars.initLowerBound(stepVars.zeroIdx:(stepVars.zeroIdx+lengthOut-1)) = lowerBound;
                    stepVars.initUpperBound(stepVars.zeroIdx:(stepVars.zeroIdx+lengthOut-1)) = upperBound;
                end
                
                if inter
                    stepVars.interLowerBound(stepVars.zeroIdx:(stepVars.zeroIdx+lengthOut-1)) = lowerBound;
                    stepVars.interUpperBound(stepVars.zeroIdx:(stepVars.zeroIdx+lengthOut-1)) = upperBound;
                end
                
                if final
                    stepVars.finalLowerBound(stepVars.zeroIdx:(stepVars.zeroIdx+lengthOut-1)) = lowerBound;
                    stepVars.finalUpperBound(stepVars.zeroIdx:(stepVars.zeroIdx+lengthOut-1)) = upperBound;
                end
                    
            case 'cost'
                stepVars.cost(stepVars.zeroIdx:(stepVars.zeroIdx+lengthOut-1)) = 1;
            case 'subSysInport'
                % points to original variable
                sameOptIndex = varargin{1};
                stepVars.optVarIdx(stepVars.zeroIdx:(stepVars.zeroIdx+lengthOut-1)) = ...
                    (sameOptIndex):(sameOptIndex+lengthOut-1);
            case 'from'
                % points to the original outport
                sameOptIndex = varargin{1};
                stepVars.optVarIdx(stepVars.zeroIdx:(stepVars.zeroIdx+lengthOut-1)) = ...
                    (sameOptIndex):(sameOptIndex+lengthOut-1);
            case 'intoSubsys'    
                % points to the original outport
                sameOptIndex = varargin{1};
                stepVars.optVarIdx(stepVars.zeroIdx:(stepVars.zeroIdx+lengthOut-1)) = ...
                    (sameOptIndex):(sameOptIndex+lengthOut-1);
            case 'demux'
                % the outport that goes into the demux points to the demux
                % outports. NOTE: This assumes that 
                sameOptIndex = varargin{1};
                outportNumber = varargin{2};
                stepVars.optVarIdx(stepVars.zeroIdx:(stepVars.zeroIdx+lengthOut-1)) =...
                    (sameOptIndex+outportNumber-1);
            case 'mux'
                % points to the outports that go into the mux
                sameOptIndex = varargin{1};
                stepVars.optVarIdx(stepVars.zeroIdx:(stepVars.zeroIdx+lengthOut-1)) = sameOptIndex;
            case 'rememberIndex'
                varargout{1} = stepVars.zeroIdx;
            case 'unitDelay'
                stepVars.state(stepVars.zeroIdx:(stepVars.zeroIdx+lengthOut-1)) = true;
            case 'input'
                stepVars.input(stepVars.zeroIdx:(stepVars.zeroIdx+lengthOut-1)) = true;
            case 'external'
                stepVars.external(stepVars.zeroIdx:(stepVars.zeroIdx+lengthOut-1)) = true;
            case 'normal'
                % do nothing
            otherwise
                % do nothing
        end
        
        % populate block.stepOutputIdx with index of the first stepVars
        % variable. Does not take into consideration redundancies. This
        % will be taken care of later.
        block.stepOutputIdx{currentBlockIndex}(portNumber) = stepVars.zeroIdx;
        
        stepVars.zeroIdx = stepVars.zeroIdx+lengthOut;
    else
        dimension = get_param(currentOutport,'CompiledPortDimensions');
        lengthOut = dimension(1)*dimension(2);
        switch stepVarsState
            case 'bound'
                fromIndex = find(stepVars.outportHandle==currentOutport,1);
                % sometimes there will be two bound blocks connected to the
                % same outport. This means that the bounds are different at
                % different time steps
                boundHandle = varargin{1};
                lowerBound = eval(get_param(boundHandle,'lb'));
                upperBound = eval(get_param(boundHandle,'ub'));

                % find out for which time steps these bounds are relevant
                init = strcmp(get_param(boundHandle, 'initial_step'), 'on');
                final = strcmp(get_param(boundHandle, 'final_step'), 'on');
                inter = strcmp(get_param(boundHandle, 'intermediate_step'), 'on');
                
                if init
                    stepVars.initLowerBound(fromIndex:(fromIndex+lengthOut-1)) = lowerBound;
                    stepVars.initUpperBound(fromIndex:(fromIndex+lengthOut-1)) = upperBound;
                end
                
                if inter
                    stepVars.interLowerBound(fromIndex:(fromIndex+lengthOut-1)) = lowerBound;
                    stepVars.interUpperBound(fromIndex:(fromIndex+lengthOut-1)) = upperBound;
                end
                
                if final
                    stepVars.finalLowerBound(fromIndex:(fromIndex+lengthOut-1)) = lowerBound;
                    stepVars.finalUpperBound(fromIndex:(fromIndex+lengthOut-1)) = upperBound;
                end
            case 'rememberIndex'
            % if this outport has already been found but we need the index for
            % the sameOptVar field
                varargout{1} = find(stepVars.outportHandle==currentOutport,1);
            case 'from'
                % points to the original outport
                fromIndex = find(stepVars.outportHandle==currentOutport,1);
                sameOptIndex = varargin{1};
                stepVars.optVarIdx(fromIndex:(fromIndex+lengthOut-1)) = ...
                    (sameOptIndex):(sameOptIndex+lengthOut-1);
            case 'subSysInport'
                % points to the original outport
                inportIndex = find(stepVars.outportHandle==currentOutport,1);
                sameOptIndex = varargin{1};
                stepVars.optVarIdx(inportIndex:(inportIndex+lengthOut-1)) = ...
                    (sameOptIndex):(sameOptIndex+lengthOut-1);
            case 'intoSubsys'    
                % points to the original outport
                subsysIndex = find(stepVars.outportHandle==currentOutport,1);
                sameOptIndex = varargin{1};
                stepVars.optVarIdx(subsysIndex:(subsysIndex+lengthOut-1)) = ...
                    (sameOptIndex):(sameOptIndex+lengthOut-1);
                
            case 'demux'
                % the outport that goes into the demux points to the demux
                % outports. NOTE: This assumes that 
                varIndex = find(stepVars.outportHandle==currentOutport,1);
                sameOptIndex = varargin{1};
                outportNumber = varargin{2};
                stepVars.optVarIdx(varIndex:(varIndex+lengthOut-1)) =...
                    (sameOptIndex+outportNumber-1);
            case 'mux'
                % points to the outports that go into the mux
                sameOptIndex = varargin{1};
                varIndex = find(stepVars.outportHandle==currentOutport,1);
                stepVars.optVarIdx(varIndex:(varIndex+lengthOut-1)) = sameOptIndex;
                
            case 'unitDelay'
                varIndex = find(stepVars.outportHandle==currentOutport,1);
                stepVars.state(varIndex:(varIndex+lengthOut-1)) = true;
                
            case 'input'
                varIndex = find(stepVars.outportHandle==currentOutport,1);
                stepVars.input(varIndex:(varIndex+lengthOut-1)) = true;  
                
            case 'external'
                varIndex = find(stepVars.outportHandle==currentOutport,1);
                stepVars.external(varIndex:(varIndex+lengthOut-1)) = true;
                
            otherwise
                % do nothing
        end
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
    referenceBlock = get_param(currentBlockHandle,'ReferenceBlock');
    
    
    % if the size of block equals the block.zeroIdx, need to double to
    % length of block
    if block.zeroIdx == length(block.handles)
        for field={'names', 'P','K','stepInputIdx','stepOutputIdx','dimensions'}
                block.(field{1}) = [block.(field{1}); cell(block.zeroIdx,1)];
        end
        block.handles = [block.handles; zeros(block.zeroIdx,1)];
        block.bound = [block.bound; false(block.zeroIdx,1)];
        block.cost = [block.cost; false(block.zeroIdx,1)];

    end
    

    if ~any(block.handles==currentBlockHandle)
        % no duplicate blocks, add this block
        currentBlockName = get_param(currentOutport,'Parent');
        block.names{block.zeroIdx} = currentBlockName;
        block.handles(block.zeroIdx) = currentBlockHandle;
        if strcmp(referenceBlock,'BLOM_Lib/Polyblock')
            % store P&K matricies if the current block is a polyblock
            block.P{block.zeroIdx} = eval(get_param(currentBlockHandle,'P'));
            block.K{block.zeroIdx}= eval(get_param(currentBlockHandle,'K'));
        elseif strcmp(referenceBlock, 'BLOM_Lib/Bound')
            block.bound(block.zeroIdx) = true;
        elseif strcmp(referenceBlock, 'BLOM_Lib/DiscreteCost')
            block.cost(block.zeroIdx) = true;
        else
            % store P and K matricies for the other blocks, as well as
            % boolean 1 if special function required and 0 ottherwise
            [P,K,specialFunPresence] = BLOM_Convert2Polyblock(currentBlockHandle);
            block.P{block.zeroIdx} = P;
            block.K{block.zeroIdx}= K;
            block.specialFunPresence{block.zeroIdx}=specialFunPresence;
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
    % first get truncated P and K matrices with the relevant times and
    % variables
    [initP,initK] = trimPK(stepP,stepK,stepVars.initTime);
    [interP,interK] = trimPK(stepP,stepK,allvars.interTime);
    [finalP,finalK] = trimPK(stepP,stepK,stepVars.finalTime);
    
    % create the P and K matrices for all intermediate time steps
    interP_full = kron(speye(horizon-2),interP);
    interK_full = kron(speye(horizon-2),interK);
    
    fullP = blkdiag(initP,interP_full,finalP);
    allK = blkdiag(initK,interK_full,finalK);
    
    % use allVars.PKOptVarIdxReroute to reroute columns of states
    % NOTE: This only works if there's only one rerouting needed. In other
    % words, only one column i can add to a certain column j. Another
    % column k cannot added to that same column j. 
    [rowPLength,~] = size(allP)
    
    % create a matrix that finds which columns of the entire P (with
    % states) to add to the columns of P without states
    toAddCol = sparse(1:rowPLength,allVars.PKOptVarIdxReroute);
    
    % create the final allP that has proper states.
    allP = toAddCol*fullP;
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
    removeP_rows = ~(any(stepP(:,~relevantTimes)~=0,2));
    trimP(removeP_rows,:) = [];
    % remove the corresponding columns of K
    trimK = stepK(:,removeP_rows);
    % remove rows of K that had a value in the removed columns of K in the
    % previous line
    removeK_rows = ~any(stepK(:,~removeP_rows)~=0,2);
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
%==========================================================================
%> @brief checks for cycles within vector (currently slower than check as you go)
%>
%> @param vector v you are trying to check
%>
%> @retval hasCycle boolean for whether a cycle exists in vector
%==========================================================================
function hasCycle = checkForCycle(v)
    hasCycle = false;
    visited = false(size(v));
    
    for i = 1:length(v)
        if (v(i) == i)
            visited(i) = true;
        end
    end

    for i = 1:length(v)
        target = i;
        traversedTargets = zeros(size(v));
        k = 1;
        while (~visited(target) || any(traversedTargets == target))
            if(any(traversedTargets == target))
                hasCycle = true;
                return
            end
            visited(target) = true;
            traversedTargets(k) = target;
            k = k + 1;
            target = v(target);
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
    
    blockIdxs = 1:block.zeroIdx-1;
    
    startBlock = [blockIdxs(block.bound | block.cost) zeros(1,10)];
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
                blockType = get_param(blockHandle, 'BlockType');
                refBlock = get_param(blockHandle, 'ReferenceBlock');
                
                if strcmp(blockType, 'SubSystem') && isempty(refBlock)
                    outportBlocks = find_system(blockHandle,'SearchDepth',1,'regexp','on','BlockType','Outport');
                    inputs = [block.stepInputIdx{blocks(idx)}];
                    for n = 1:length(outportBlocks)
                        outportBlockPorts = get_param(outportBlocks(n), 'PortHandles');
                        outportInputHandle = outportBlockPorts.Inport;
                        line = get_param(outportInputHandle, 'Line');
                        srcPortHandle = get_param(line, 'SrcPortHandle');
                        inputsToAdd = find(stepVars.outportHandle == srcPortHandle);
                        inputs = [inputs; inputsToAdd];
                    end

                elseif strcmp(blockType, 'From')
                    gotoTag = get_param(blockHandle, 'GotoTag');
                    name = get_param(blockHandle, 'Parent');
                    gotoBlock = find_system(name,'BlockType','Goto','GotoTag',gotoTag);
                    gotoPorts = get_param(gotoBlock{1},'PortHandles');
                    inputHandle = gotoPorts.Inport;
                    line = get_param(inputHandle, 'Line');
                    srcPortHandle = get_param(line, 'SrcPortHandle');
                    inputs = find(stepVars.outportHandle == srcPortHandle);
                    
                elseif strcmp(blockType, 'UnitDelay') && (idx ~= 1)
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
    
    % label what timestep each variable exists in
    allVars.timeStep(1:initialLength) = ones(initialLength,1);
    interSteps = ones(interLength,1) * (2:horizon-1);
    allVars.timeStep(initialLength+1:end-finalLength) = interSteps(:);
    allVars.timeStep(end-finalLength+1:end) = horizon * ones(finalLength,1);
    
    
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
        kron(ones(horizon-2,1),stepVarsIndices(stepVars.interTime));
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
    
    OptVarIdxReroute = zeros(max(allVars.optVarIdx),1);
    allVars.PKOptVarIdxReroute = zeros(max(allVars.optVarIdx),1);
    
    for idx = initialLength+1:totalLength
        stepVarIdx = allVars.stepVarIdx(idx);
        if (stepVars.state(stepVarIdx))
            blockIdx = stepVars.block(stepVarIdx);
            timeStep = allVars.timeStep(idx);
            blockInputOutputIdx = find(block.stepOutputIdx{blockIdx}== stepVarIdx, 1);
            newOptVarIdx = allVars.optVarIdx(block.allInputMatrix{blockIdx}(blockInputOutputIdx,timeStep-1));
            
            OptVarIdxReroute(allVars.optVarIdx(idx)) = newOptVarIdx;
            allVars.optVarIdx(idx) = newOptVarIdx;
        end
    end

    selfMapIdxs = PKOptVarIdxReroute == 0;
    allVars.PKOptVarIdxReroute(selfMapIdxs) = 1:length(selfMapIdxs);
    for stateMapIdx = find(PKOptVarIdxReroute ~= 0)
        allVars.PKOptVarIdxReroute(stateMapIdx) = allVars.PKOptVarIdxReroute(OptVarIdxReroute(stateMapIdx));
    end

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
function [ModelSpec] = convert2ModelSpec(ModelSpec,stepVars,block,varargin)
    ModelSpec.in_vars = stepVars.input;
    ModelSpec.ex_vars = stepVars.external;
end