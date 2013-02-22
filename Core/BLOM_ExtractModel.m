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

function [ModelSpec,block,allVars] = BLOM_ExtractModel(name,horizon,dt,integ_method,options)
    % load system. does nothing if model is not open
    load_system(name);
    % evaluate model to get dimensions
    eval([name '([],[],[],''compile'');']); 
    try
        [boundHandles,costHandles,inputAndExternalHandles] = findBlocks(name);
        [block,allVars,stop] = ...
            searchSources(boundHandles,costHandles,inputAndExternalHandles,name);
        % FIX: should implement something that stops the code after analyzing
        % all the blocks and finding an error in the structure of the model
        if stop == 1
            % break the code somehow?
        end

        % find out which wires are relevant at which times
        %[timeStruct] = relevantTimes(outportHandles);
        
        % create large P and K matrix
       % try
          %  [bigP,bigK] = combinePK(block,allVars);
        %catch err
           % rethrow(err)
        %end
        
        % following code checks whether or not inports and outportHandles
        % was filled in properly
        fprintf('--------------------------------------------------------\n')
        fprintf('Test block.outputIdxs and inports is correct\n')
        fprintf('Everything is good if nothing prints\n')
        fprintf('--------------------------------------------------------\n')
        for i = 1:length(block.handles)
            currentPorts = get_param(block.handles(i),'PortHandles');
            currentOutports = currentPorts.Outport;
            currentInports = currentPorts.Inport;
            % compare lists of Outports
            
            currentBlockOutports = zeros(length(block.outputIdxs{i}),1);
            for idx = 1:length(block.outputIdxs{i})
                currentBlockOutports(idx) = allVars.outportHandle(block.outputIdxs{i}(idx));
            end
            
            diffOutports = setdiff(currentBlockOutports,currentOutports);
            if ~isempty(diffOutports)
                fprintf('Difference in Outports in %s\n',block.names{i})
                currentOutports
                get_param(allVars.outportHandle(block.outputIdxs{i}),'Parent')
                fprintf('--------------------------------------------------------\n')
            end
            
            % compare inports
            inportsOutport = zeros(length(currentInports),1);
            for index = 1:length(currentInports)
                currentLine = get_param(currentInports(index),'Line');
                inportsOutport(index) = get_param(currentLine,'SrcPortHandle');
            end
            
            currentBlockInports = zeros(length(block.inputIdxs{i}),1);
            for idx = 1:length(block.inputIdxs{i})
                currentBlockInports(idx) = allVars.outportHandle(block.inputIdxs{i}(idx));
            end
            
            diffInports = setdiff(inportsOutport,currentBlockInports);
            if ~isempty(diffInports) && ~any(block.handles(i)==inputAndExternalHandles)
                fprintf('Difference in inputs in %s\n',block.names{i})
                inportsOutport
                block.inputIdxs{i}
                get_param(inportsOutport,'Parent')
                fprintf('--------------------------------------------------------\n')
            end
            
            
            % first need to get list of inports for each block and the
            % outport that's associated with that inport
            
        end
        
        %following code is to make sure searchSources works
        fprintf('The Number of blocks found is %.0f\n',length(block.handles))
        fprintf('The Number of outports found is %.0f\n',length(allVars.outportHandle))
        for i = 1:length(allVars.outportHandle);
            parent = get_param(allVars.outportHandle(i),'Parent');
            allVars.outportHandle(i);
            portType = get_param(allVars.outportHandle(i),'PortType');
            if ~strcmp(portType,'outport')
                fprintf('Oops, this is not an outport, look to see what happened\n')
                portType
            end
        end
        

        fprintf('\n\n\n Here we have the stored data from block\n\n\n');

        for i = 1:length(block.handles);
            someBlock = block.names{i};
        end
        
        % check to see if allVars points to the proper block
        fprintf('Test to make sure allVars points to the proper block\n')
        fprintf('--------------------------------------------------------\n')
        for i = 1:length(allVars.block)
            currentBlockName = block.names(allVars.block(i));
            parent = get_param(allVars.outportHandle(i),'Parent');
            
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
        [...find_system(name,'FindAll','on','ReferenceBlock','BLOM_Lib/InputFromWorkspace'); ...
        find_system(name,'FindAll','on','ReferenceBlock','BLOM_Lib/InputFromSimulink'); ...
        ...find_system(name,'FindAll','on','ReferenceBlock','BLOM_Lib/ExternalFromWorkspace'); ...
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
%> @retval allVars structure that contains information about each
%> optimization variable
%> @retval stop if there are any unconnected blocks or blocks that BLOM
%> does not support, this parameter gives a 1 to indicate true and 0 to
%> indicate false
%======================================================================

function [block,allVars,stop] = searchSources(boundHandles,costHandles,...
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
    
    allVarsZero = 1; % index of first allVars zero
    % allVars stores information for every optimization variable
    allVars.block = zeros(initialSize,1); % the index of the block
    allVars.outportNum = zeros(initialSize,1); % outport number. e.g. if there are multiple outports of a block
    allVars.outportHandle = zeros(initialSize,1); % handle of specific outport
    allVars.outportIndex = zeros(initialSize,1); % index of specific outport. normally just 1 (will be more than 1 if the outport is a vector)
    allVars.optVarIdx = (1:initialSize)'; % points to the true optimization variable index. (will usually point to itself, otherwise, some "true" variable)
    allVars.cost = zeros(initialSize,1); %default cost is zero. 1 if we see that it's part of the cost
    allVars.upperBound = inf*ones(initialSize,1); % upper bound of each variable. default inf
    allVars.lowerBound = -inf*ones(initialSize,1); % lower bound of each variable. default -inf
    allVars.time = cell(initialSize,1);
    
    blockZero = 1; % index of first block zero
    % block stores information about each block
    block.names = cell(initialSize,1); % name of block
    block.handles = zeros(initialSize,1); % handle of block
    block.P = cell(initialSize,1); % P matrix of block, if relevant
    block.K = cell(initialSize,1); % K matrix of block, if relevant
    block.inputIdxs = cell(initialSize,1); % the first optVarIDx of the outports connected to the inports of that block
    block.outputIdxs = cell(initialSize,1); %the first optVarIdx of the outports of that block
    block.dimensions = cell(initialSize,1); % dimensions of each outport. first value is outport #, then second two values are dimensions of outport
    block.sourceOutports = cell(initialSize,1); % source outport handles
    
    % find all lines connected to costs and bounds and then get outport
    % ports from there. fill in allVars and block structures as necessary
    
    % iZero is index of first zero of outportHandles
    iZero = 1;
    outportHandles = zeros(initialSize,1);
    
    % get all outports connected to the bounds 
    state = 'bound';
    for i = 1:length(boundHandles)
        portH = get_param(boundHandles(i),'PortHandles');
        % costs and bounds should only have one inport and line
        currentInport = [portH.Inport];
        [outportHandles,iZero,allVars,allVarsZero,block,blockZero] = ...
            updateVars(currentInport,outportHandles,iZero,...
            allVars,allVarsZero,block,blockZero,state,iZero,portH,boundHandles(i));
    end
 
    
    % get all outports connected to costs
    state = 'cost';
    for i = 1:length(costHandles)
        portH = get_param(costHandles(i),'PortHandles');
        % costs and bounds should only have one inport and line
        currentInports = [portH.Inport];
        [outportHandles,iZero,allVars,allVarsZero,block,blockZero] = ...
            updateVars(currentInports,outportHandles,iZero,...
            allVars,allVarsZero,block,blockZero,state);
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
            [outportHandles,iZero,allVars,allVarsZero,block,blockZero] = ...
                    updateVars(sourceInports,outportHandles,iZero,...
                    allVars,allVarsZero,block,blockZero,state,iOut,sourcePorts);
            
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
                [outportHandles,iZero,allVars,allVarsZero,block,blockZero] = ...
                    updateVars(sourceOutports,outportHandles,iZero,...
                    allVars,allVarsZero,block,blockZero,state,iOut,sourcePorts);

            elseif strncmp(refBlock,'BLOM_Lib',8)
                % subsystem that's a BLOM block
                state = 'BLOM';
                [outportHandles,iZero,allVars,allVarsZero,block,blockZero] = ...
                    updateVars(sourceInports,outportHandles,iZero,...
                    allVars,allVarsZero,block,blockZero,state,iOut,sourcePorts);
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
            [outportHandles,iZero,allVars,allVarsZero,block,blockZero] = ...
                updateVars(sourceInports,outportHandles,iZero,...
                allVars,allVarsZero,block,blockZero,state,iOut,sourcePorts);
            
        elseif strcmp(sourceType,'Demux')
            % block is a demux. Here we want to fill in field sameOptVar to
            % point to the original variable
            
            state = 'demux';
            sourceOutports = [sourcePorts.Outport];
            [outportHandles,iZero,allVars,allVarsZero,block,blockZero] = ...
                updateVars(sourceInports,outportHandles,iZero,...
                allVars,allVarsZero,block,blockZero,state,iOut,sourcePorts,sourceOutports);
            
        elseif strcmp(sourceType,'Mux')
            % block is a mux. Here we want to fill in field sameOptVar to
            % point to the original variable
            
            state='mux';
            [outportHandles,iZero,allVars,allVarsZero,block,blockZero] = ...
                updateVars(sourceInports,outportHandles,iZero,...
                allVars,allVarsZero,block,blockZero,state,iOut,sourcePorts);
            
        elseif length(sourceInports) > 1
            % currently do nothing special if there is more than one input
            % except look for what's connected. 
            % FIX: check to see which inports affect which outports. Not
            % all inports may be relevant
            state = 'norm';
            [outportHandles,iZero,allVars,allVarsZero,block,blockZero] = ...
                updateVars(sourceInports,outportHandles,iZero,...
                allVars,allVarsZero,block,blockZero,state,iOut,sourcePorts);
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
                        [outportHandles,iZero,allVars,allVarsZero,block,blockZero] = ...
                            updateVars(sourceInports(portNumber),outportHandles,iZero,...
                            allVars,allVarsZero,block,blockZero,state,iOut,sourcePorts);
                    else
                        % This case should NOT be reached because of an
                        % inport exists, there is almost certainly a
                        % corresponding inport for the entire subsystem
                        % itself. It is here as asafeguard. 
                        allVarsState = 'normal';
                        fprintf('Have you forgotten to connect one of the inports to subsystem: %s? ', parentOfBlock);
                    end
                else
                    allVarsState = 'normal';
                    [block,blockZero,currentBlockIndex] =...
                        updateBlock(block,blockZero,outportHandles(iOut));
                    [allVars,allVarsZero,block] = updateAllVars(allVars,allVarsZero,...
                        block,currentBlockIndex,outportHandles(iOut),allVarsState);
                end
            else
                allVarsState = 'normal';
                [block,blockZero,currentBlockIndex] =...
                    updateBlock(block,blockZero,outportHandles(iOut));
                [allVars,allVarsZero,block] = updateAllVars(allVars,allVarsZero,...
                    block,currentBlockIndex,outportHandles(iOut),allVarsState);
                iOut = iOut+1;
                continue
            end
        else
            % in the case of one inport, it must affect the outport, so
            % find the relevant outports
            [outportHandles,iZero,allVars,allVarsZero,block,blockZero] = ...
                updateVars(sourceInports,outportHandles,iZero,...
                allVars,allVarsZero,block,blockZero,2,iOut,sourcePorts);
        end
        % FIX: should check to see if BLOM supports the blocks that is
        % found right here
        iOut = iOut+1;
    end
    
    
    if any(outportHandles==0)
        outportHandles = outportHandles(1:iZero-1);
    end
    
    block = updateInputsField(block,allVars,outportHandles);
    
    %check if any handles are -1 and remove them
    if any(outportHandles==-1)
        % FIX PRINT STATEMENT. NOT NECESSARILY TRUE
        fprintf('One or more of the blocks are missing a connection\n');
        stop = 1;
    end
    
    % remove empty and 0 entries in blocks
    for field={'names', 'P','K','inputIdxs','outputIdxs','dimensions',...
            'sourceOutports'}
        block.(field{1}) = block.(field{1})(1:(blockZero-1));
    end
    block.handles = block.handles(1:(blockZero-1));
    
    % remove empty and 0 entries in allVars
    for field={'block', 'outportNum','outportHandle','outportIndex',...
        'optVarIdx','cost'}
        allVars.(field{1}) = allVars.(field{1})(1:(allVarsZero-1));
    end
    allVars.upperBound = allVars.upperBound(1:(allVarsZero-1));
    allVars.lowerBound = allVars.lowerBound(1:(allVarsZero-1));
    allVars.time = allVars.time(1:(allVarsZero-1));
    % FIX: need to find some way to remove all -1 handles. using setdiff
    % with [-1] reorders all the outport handles and puts it in ascending
    % order
    % outportHandles = setdiff(outportHandles,[-1]);
    
end

%======================================================================
%> @brief From given inports, see which outports are relevant. Also, given
%> the current outport, populate the block and allVars structures
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

function [outportHandles,iZero,allVars,allVarsZero,block,blockZero] =...
    updateVars(inports,existingOutports,iZero,allVars,allVarsZero,...
    block,blockZero,state,varargin)

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
            allVarsState = 'bound';
            boundHandle = varargin{3};
            lowerBound = eval(get_param(boundHandle,'lb'));
            upperBound = eval(get_param(boundHandle,'ub'));
            
            % go through each of the outports connected through bound.
            % Because bound has only one inport, the for loop above only
            % goes through one iteration, so the variable outportFound
            % will include all the outports found
            for currentOutport = outportFound
                [block,blockZero,currentBlockIndex] =...
                    updateBlock(block,blockZero,currentOutport);
                [allVars,allVarsZero,block] = updateAllVars(allVars,allVarsZero,...
                    block,currentBlockIndex,currentOutport,allVarsState,lowerBound,upperBound);
            end
            
        case 'cost'
            allVarsState = 'cost';
            % go through each of the outports connected through cost
            % Because cost has only one inport, the for loop above only
            % goes through one iteration, so the variable outportFound
            % will include all the outports found
            for currentOutport = outportFound
                [block,blockZero,currentBlockIndex] =...
                    updateBlock(block,blockZero,currentOutport);
                [allVars,allVarsZero,block] = updateAllVars(allVars,allVarsZero,...
                    block,currentBlockIndex,currentOutport,allVarsState);
            end
            
        case 'fromSimulink'
            % for the fromSimulink blocks, we need a special case to look
            % at that specific outport even though we stop the BFS here
            allVarsState = 'normal';
            currentOutport = existingOutports(iOut);
            [block,blockZero,currentBlockIndex] =...
                updateBlock(block,blockZero,currentOutport);
            [allVars,allVarsZero,block] = updateAllVars(allVars,allVarsZero,...
                block,currentBlockIndex,currentOutport,allVarsState);
            
        case 'subSysInport'
            % Here we fill in the optVarIdx field have the inport point to
            % the outport connected to the inport of the subsystem

            % the current outport changes here because we want the outport
            % connected to the subsystem which is connected to the inport
            % to the subsystem to have a reference that it's a duplicate
            % variable
            allVarsState = 'rememberIndex';
            currentOutport = outportFound;
            [block,blockZero,currentBlockIndex] =...
                updateBlock(block,blockZero,currentOutport);
            [allVars,allVarsZero,block,sameOptIndex] = updateAllVars(allVars,allVarsZero,...
                block,currentBlockIndex,currentOutport,allVarsState);
            
            allVarsState = 'subSysInport';
            currentOutport = existingOutports(iOut);
            [block,blockZero,currentBlockIndex] =...
                updateBlock(block,blockZero,currentOutport);
            [allVars,allVarsZero,block] = updateAllVars(allVars,allVarsZero,...
                    block,currentBlockIndex,currentOutport,allVarsState,sameOptIndex);
        case 'from'
            % in this case, we want to be able to fill in the sameOptVar
            % field in allVars so that we don't create a separate
            % optimization variable but rather know that this is a duplicate
            allVarsState = 'rememberIndex';
            % the currentOutport in this case is the outport that goes to
            % the GoTo block
            currentOutport = outportFound;
            [block,blockZero,currentBlockIndex] =...
                updateBlock(block,blockZero,currentOutport);
            [allVars,allVarsZero,block,sameOptIndex] = updateAllVars(allVars,allVarsZero,...
                    block,currentBlockIndex,currentOutport,allVarsState);
            % save information for from block here
            allVarsState = 'from';
            currentOutport = existingOutports(iOut);
            [block,blockZero,currentBlockIndex] =...
                updateBlock(block,blockZero,currentOutport);
            [allVars,allVarsZero,block] = updateAllVars(allVars,allVarsZero,...
                    block,currentBlockIndex,currentOutport,allVarsState,sameOptIndex);
                
        case 'intoSubsys'
            % for this case, we want to fill in the sameOptVar variable for
            % the subsystems outport and say the original variable is the
            % outport of whatever is inside the subsystem
            allVarsState = 'rememberIndex';
            % the currentOutport in this case is the outport that goes to
            % the subsystem
            currentOutport = outportFound;
            [block,blockZero,currentBlockIndex] =...
                updateBlock(block,blockZero,currentOutport);
            [allVars,allVarsZero,block,sameOptIndex] = updateAllVars(allVars,allVarsZero,...
                block,currentBlockIndex,currentOutport,allVarsState);
            
            % save information for subsystem outport here
            allVarsState = 'intoSubsys';
            currentOutport = existingOutports(iOut);
            [block,blockZero,currentBlockIndex] =...
                updateBlock(block,blockZero,currentOutport);
            [allVars,allVarsZero,block] = updateAllVars(allVars,allVarsZero,...
                block,currentBlockIndex,currentOutport,allVarsState,sameOptIndex);
                
        case 'demux'
            % the original variable is going to be the outport connected to
            % the demux
            
            allVarsState = 'rememberIndex';
            % should only be one outport connected to the demux.
            currentOutport = outportFound;
            [block,blockZero,currentBlockIndex] =...
                updateBlock(block,blockZero,currentOutport);
            [allVars,allVarsZero,block,sameOptIndex] = updateAllVars(allVars,allVarsZero,...
                block,currentBlockIndex,currentOutport,allVarsState);
            
            % fill in optVarIdx for the demux's outports to the outport
            % that goes into demux
            sourceOutports = varargin{3};
            for i = 1:length(sourceOutports)
                % store information for the rest of the demux outports
                allVarsState = 'demux';
                [block,blockZero,currentBlockIndex] =...
                    updateBlock(block,blockZero,sourceOutports(i));
                [allVars,allVarsZero,block] = updateAllVars(allVars,allVarsZero,...
                    block,currentBlockIndex,sourceOutports(i),allVarsState,sameOptIndex,i);
            end
                
            
            
        case 'mux'           
            % the original variables are the outports connected to the mux
            mux_optVarIdx = zeros(length(allOutportsFound),1);
            for idx=1:length(allOutportsFound)
                   currentOutport = allOutportsFound(idx);                
                   allVarsState = 'rememberIndex';
                    [block,blockZero,currentBlockIndex] =...
                        updateBlock(block,blockZero,currentOutport);
                    [allVars,allVarsZero,block,sameOptIndex] = updateAllVars(allVars,allVarsZero,...
                        block,currentBlockIndex,currentOutport,allVarsState);
                    mux_optVarIdx(idx) = sameOptIndex;
             end    
              % save information for mux outport here   
                allVarsState = 'mux';
                currentOutport=existingOutports(iOut);
                [block,blockZero,currentBlockIndex] =...
                    updateBlock(block,blockZero,currentOutport);
                [allVars,allVarsZero,block] = updateAllVars(allVars,allVarsZero,...
                    block,currentBlockIndex,currentOutport,allVarsState,mux_optVarIdx);  

        otherwise 
            % not looking at bounds or costs
            allVarsState = 'normal';
            %update block structure
            [block,blockZero,currentBlockIndex] =...
                updateBlock(block,blockZero,currentOutport);
            [allVars,allVarsZero,block] = updateAllVars(allVars,allVarsZero,...
                block,currentBlockIndex,currentOutport,allVarsState);
    end
        
    
end

%%
%======================================================================
%> @brief update allVars structure 
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
%> @param allVars allVars structure
%> @param allVarsZero current index of first zero of allVars
%> @param blockZero current index of first zero of block structure. used to
%> figure out which block the outport is from
%> @param currentOutport outport that you want to add information for
%> @param state lets the function know how to populate allVars. current
%> states are 'bound', 'cost'
%>
%> @retval allVars allVars structure
%> @retval allVarsZero updated index of first zero of allVarsZero
%======================================================================

function [allVars,allVarsZero,block,varargout] = updateAllVars(allVars,allVarsZero,...
    block,currentBlockIndex,currentOutport,allVarsState,varargin)
%this function populates allVars structure
%state can be 'bound', 'cost' 
    if sum(any(allVars.outportHandle==currentOutport)) == 0
        %only add outports which haven't been looked at
        dimension = get_param(currentOutport,'CompiledPortDimensions');
        lengthOut = dimension(1)*dimension(2);
        currentIndices = 1:lengthOut;
        portNumber = get_param(currentOutport,'PortNumber');
        
        % update the size of allVars as necessary
        newLength = allVarsZero+lengthOut;
        oldLength = length(allVars.block);
        if newLength >= oldLength
            if newLength >= 2*oldLength
            % new entries greater than twice the old length
                for field={'block', 'outportNum','outportHandle','outportIndex',...
                        'cost'}
                        allVars.(field{1}) = [allVars.(field{1}); zeros(newLength*2,1)];
                end
                allVars.optVarIdx = [allVars.optVarIdx; ((oldLength+1):(newLength*2+oldLength))'];
                allVars.upperBound = [allVars.upperBound; inf*ones(newLength*2,1)];
                allVars.lowerBound = [allVars.lowerBound; -inf*ones(newLength*2,1)];
                allVars.time = [allVars.time; cell(newLength*2,1)];
            else % double the length of all fields in allVars
                for field={'block', 'outportNum','outportHandle','outportIndex',...
                        'cost'}
                        allVars.(field{1}) = [allVars.(field{1}); zeros(oldLength,1)];
                end
                allVars.optVarIdx = [allVars.optVarIdx; ((oldLength+1):(oldLength*2))'];
                allVars.upperBound = [allVars.upperBound; inf*ones(oldLength,1)];
                allVars.lowerBound = [allVars.lowerBound; -inf*ones(oldLength,1)];
                allVars.time = [allVars.time; cell(oldLength,1)];
            end
        end
        allVars.block(allVarsZero:(allVarsZero+lengthOut-1)) = currentBlockIndex;
        allVars.outportNum(allVarsZero:(allVarsZero+lengthOut-1)) = portNumber;
        allVars.outportHandle(allVarsZero:(allVarsZero+lengthOut-1)) = currentOutport;
        allVars.outportIndex(allVarsZero:(allVarsZero+lengthOut-1)) = currentIndices;
        
        switch allVarsState
            case 'bound'
                lowerBound = varargin{1};
                upperBound = varargin{2};
                allVars.lowerBound(allVarsZero:(allVarsZero+lengthOut-1)) = lowerBound;
                allVars.upperBound(allVarsZero:(allVarsZero+lengthOut-1)) = upperBound;
            case 'cost'
                allVars.cost(allVarsZero:(allVarsZero+lengthOut-1)) = 1;
            case 'subSysInport'
                % points to original variable
                sameOptIndex = varargin{1};
                allVars.optVarIdx(allVarsZero:(allVarsZero+lengthOut-1)) = ...
                    (sameOptIndex):(sameOptIndex+lengthOut-1);
            case 'from'
                % points to the original outport
                sameOptIndex = varargin{1};
                allVars.optVarIdx(allVarsZero:(allVarsZero+lengthOut-1)) = ...
                    (sameOptIndex):(sameOptIndex+lengthOut-1);
            case 'intoSubsys'    
                % points to the original outport
                sameOptIndex = varargin{1};
                allVars.optVarIdx(allVarsZero:(allVarsZero+lengthOut-1)) = ...
                    (sameOptIndex):(sameOptIndex+lengthOut-1);
            case 'demux'
                % the outport that goes into the demux points to the demux
                % outports. NOTE: This assumes that 
                sameOptIndex = varargin{1};
                outportNumber = varargin{2};
                allVars.optVarIdx(allVarsZero:(allVarsZero+lengthOut-1)) =...
                    (sameOptIndex+outportNumber-1);
            case 'mux'
                % points to the outports that go into the mux
                sameOptIndex = varargin{1};
                allVars.optVarIdx(allVarsZero:(allVarsZero+lengthOut-1)) = sameOptIndex;
            case 'rememberIndex'
                varargout{1} = allVarsZero;
            case 'normal'
                % do nothing
            otherwise
                % do nothing
        end
        
        % populate block.outputIdxs with index of the first allVars
        % variable. Does not take into consideration redundancies. This
        % will be taken care of later.
        block.outputIdxs{currentBlockIndex}(portNumber) = allVarsZero;
        
        allVarsZero = allVarsZero+lengthOut;
    else
        dimension = get_param(currentOutport,'CompiledPortDimensions');
        lengthOut = dimension(1)*dimension(2);
        switch allVarsState
            case 'rememberIndex'
            % if this outport has already been found but we need the index for
            % the sameOptVar field
                varargout{1} = find(allVars.outportHandle==currentOutport,1);
            case 'from'
                % points to the original outport
                fromIndex = find(allVars.outportHandle==currentOutport,1);
                sameOptIndex = varargin{1};
                allVars.optVarIdx(fromIndex:(fromIndex+lengthOut-1)) = ...
                    (sameOptIndex):(sameOptIndex+lengthOut-1);
            case 'subSysInport'
                % points to the original outport
                inportIndex = find(allVars.outportHandle==currentOutport,1);
                sameOptIndex = varargin{1};
                allVars.optVarIdx(inportIndex:(inportIndex+lengthOut-1)) = ...
                    (sameOptIndex):(sameOptIndex+lengthOut-1);
            case 'intoSubsys'    
                % points to the original outport
                subsysIndex = find(allVars.outportHandle==currentOutport,1);
                sameOptIndex = varargin{1};
                allVars.optVarIdx(subsysIndex:(subsysIndex+lengthOut-1)) = ...
                    (sameOptIndex):(sameOptIndex+lengthOut-1);
                
            case 'demux'
                % the outport that goes into the demux points to the demux
                % outports. NOTE: This assumes that 
                varIndex = find(allVars.outportHandle==currentOutport,1);
                sameOptIndex = varargin{1};
                outportNumber = varargin{2};
                allVars.optVarIdx(varIndex:(varIndex+lengthOut-1)) =...
                    (sameOptIndex+outportNumber-1);
            case 'mux'
                % points to the outports that go into the mux
                sameOptIndex = varargin{1};
                varIndex = find(allVars.outportHandle==currentOutport,1);
                allVars.optVarIdx(varIndex:(varIndex+lengthOut-1)) = sameOptIndex;
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
%> @param blockZero current index of first zero of block
%> @param currentOutport outport that you want to add information for
%>
%> @retval block block structure
%> @retval blockZero updated index of first zero of block
%======================================================================

function [block,blockZero,currentBlockIndex] = updateBlock(block,blockZero,currentOutport)

    currentBlockHandle = get_param(currentOutport,'ParentHandle');
    sourcePorts = get_param(currentBlockHandle,'PortHandles');
    referenceBlock = get_param(currentBlockHandle,'ReferenceBlock');
    
    
    % if the size of block equals the blockZero, need to double to
    % length of block
    if blockZero == length(block.handles)
        for field={'names', 'P','K','inputIdxs','outputIdxs','dimensions',...
                'sourceOutports'}
                block.(field{1}) = [block.(field{1}); cell(blockZero,1)];
        end
        block.handles = [block.handles; zeros(blockZero,1)];
    end
    

    if ~any(block.handles==currentBlockHandle)
        % no duplicate blocks, add this block
        currentBlockName = get_param(currentOutport,'Parent');
        block.names{blockZero} = currentBlockName;
        block.handles(blockZero) = currentBlockHandle;
        if strcmp(referenceBlock,'BLOM_Lib/Polyblock')
            % store P&K matricies if the current block is a polyblock
            block.P{blockZero} = eval(get_param(currentBlockHandle,'P'));
            block.K{blockZero}= eval(get_param(currentBlockHandle,'K'));
        else
            % store P and K matricies for the other blocks
            [P,K] = BLOM_Convert2Polyblock(currentBlockHandle);
            block.P{blockZero} = P;
            block.K{blockZero}= K;
        end

        if isempty(block.outputIdxs{blockZero})
            block.outputIdxs{blockZero} = zeros(length(sourcePorts.Outport),1);
        end
        

        % increase the index of block by one and populate currentBlockIndex
        currentBlockIndex = blockZero;
        blockZero = blockZero+1;
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

function [block] = updateInputsField(block,allVars,outportHandles)
    % get information for the block and port that this outport goes to
    for index = 1:length(outportHandles)
        currentOutport = outportHandles(index);
        % find the first index that corresponds to this outport in allVars
        allVarsIdx = find(allVars.outportHandle==currentOutport,1);
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
                if isempty(block.inputIdxs{inportBlockIndex})
                    block.inputIdxs{inportBlockIndex} = zeros(length(inportPorts),1);
                end
                block.inputIdxs{inportBlockIndex}(inportIndex) = allVarsIdx;
            end
        end
    end

end


%%
%======================================================================
%> @brief given block and allVars, create a large P and K matrix with the
%> indicies of allVars as the columns
%>
%> More detailed description of the problem.
%>
%> @param block block structure
%> @param allVars allVars structure
%>
%> @retval bigP P matrix of all variables
%> @retval bigK K matrix of all variables
%======================================================================

function [bigP,bigK] = combinePK(block,allVars)
    pRowSize = 30;
    kRowSize = 30;
    pRowZero = 1;
    kRowZero = 1;
    bigP = zeros(pRowSize,length(allVars.block));
    bigK = zeros(kRowSize,length(allVars.block));
    
    for idx = 1:length(block.handles)
        % go through all the blocks and fill in what the P and K matrix
        % should be
        % idx = block index

        if ~isempty(block.P{idx})
            % only start filling in the P and K matrices if there are P and
            % K matrices to fill. not all blocks will have this field
            
            % find the current P and K matrix
            currentP = block.P{idx};
            pRowLength = size(currentP,1);
            currentK = block.K{idx};
            kRowLength = size(currentK,1);
            
            % find all the relevant allVars indicies that match the
            % outports of block.inputIdxs
            inputsIndices = zeros(length(block.inputIdxs{idx}),length(allVars.block));
            for inputsIdx = 1:length(block.inputIdxs{idx})
                inputsIndices(inputsIdx,:) =...
                    allVars.outportHandle==block.inputIdxs{idx}(inputsIdx);
                if sum(inputsIndices(inputsIdx,:)) == 0
                    % This case should never happen. DELETE AFTER TESTING.
                    parentName = get_param(block.inputIdxs{idx}(inputsIdx),'Parent');
                    block.names{idx};
                    fprintf('For some reason, cannot find input for %s\n',parentName)
                end
            end
            inputsIndices = logical(inputsIndices);
            
            % find all the relevant allVars indicies that match the
            % outports of block.inputIdxs
            outputsIndices = zeros(length(block.outputIdxs),length(allVars.block));
            for outputsIdx = 1:length(block.outputIdxs{idx});
                outputsIndices(outputsIdx,:) =...
                    allVars.outportHandle==block.outputIdxs{idx}(outputsIdx);
            end
            outputsIndices = logical(outputsIndices);
            
            % FIX: Increase the size of bigP and bigK when needed
            
            % here we put the relevant parts of P and K in the right places
            currentZeroCol = 1;
            sumInput = sum(inputsIndices,2);

            for pIdx = 1:size(inputsIndices,1)
                dimOutport = sumInput(pIdx);
                full(currentP)
                full(currentP(:,currentZeroCol:(currentZeroCol+dimOutport-1)))
                size(bigP(pRowZero:(pRowZero+pRowLength-1),inputsIndices(pIdx,:)))
                size(currentP(:,currentZeroCol:(currentZeroCol+dimOutport-1)))
                bigP(pRowZero:(pRowZero+pRowLength-1),inputsIndices(pIdx,:))...
                    = currentP(:,currentZeroCol:(currentZeroCol+dimOutport-1));
                currentZeroCol = currentZeroCol + dimOutport;
            end
            
            
            pRowZero = pRowZero + pRowLength;
            kRowZero = kRowZero + kRowLength;
            
        end
    end

end

%%
%==========================================================================
%> @brief Reroutes pointers such that if System A points to System B, and 
%> System B points to System C, System A points to System C instead.
%> optVarIdx then reduced to remove gaps. e.g. [1;3;4;5]->[1;2;3;4]
%> 
%> @param optVarIdx vector to be rerouted and reduced
%> 
%> @retval optVarIdx rerouted and reduced
%==========================================================================
function optVarIdx = cleanupOptVarIdx(optVarIdx)

    if (checkForCycle(optVarIdx))
        %Error, optVarIdx contains cycle
    end

    for i = 1:length(optVarIdx)
        target = i;
        while (target ~= optVarIdx(target))
            target = optVarIdx(target);
        end
        optVarIdx(i) = target;
    end
    [~,~,optVarIdx] = unique(optVarIdx);
end

%%
%==========================================================================
%> @brief checks for cycles within vector
%>
%> @param vector you are trying to check
%>
%> @retval hasCycle boolean for whether a cycle exists in vector
%==========================================================================
function hasCycle = checkForCycle(vector)

end