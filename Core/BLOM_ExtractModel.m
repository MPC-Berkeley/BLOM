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

function [ModelSpec] = BLOM_ExtractModel(name,horizon,dt,integ_method,options)
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

        %following code is to make sure searchSources works
        fprintf('The Number of blocks is %.0f\n',length(block.handles))
        fprintf('The Number of outports is %.0f\n',length(allVars.outportHandle))
        for i = 1:length(allVars.outportHandle);
            allVars.outportHandle(i);
            parent = get_param(allVars.outportHandle(i),'Parent');
            portType = get_param(allVars.outportHandle(i),'PortType');
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
                fprintf('correct index\n')
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
        % close evaluation of models
    catch err
        rethrow(err)
    end
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
%> More detailed description of the problem.
%>
%> @param handleArray the array of handles that you want to search. These
%> are the block handles of the costs and constraints
%> @param varargin the external and input from simulink handles. the
%> sources of these blocks are not relevant to the optimization problem
%>
%> @retval outportHandles returns an array of the outport handles found
%> from BFS search
%> @retval stop if there are any unconnected blocks or blocks that BLOM
%> does not support, this parameter gives a 1 to indicate true and 0 to
%> indicate false
%> @retval boundsStruct structure with fields: 1) outport handles found, 2)
%> boolean true or false if it is connected to a bounds block
%======================================================================

function [block,allVars,stop] = searchSources(boundHandles,costHandles,varargin)
    % only flag stop = 1 when there is a bad block
    stop = 0;
    if ~isempty(varargin)
        % gets blocks of final stopping places
        fromSimulink = varargin{1};
        name = varargin{2};
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

    end
    
    initialSize = 20;
    
    allVarsZero = 1; % index of first allVars zero
    % allVars stores information for every optimization variable
    allVars.block = zeros(initialSize,1); % the index of the block
    allVars.outportNum = zeros(initialSize,1); % outport number. e.g. if there are multiple outports of a block
    allVars.outportHandle = zeros(initialSize,1); % handle of specific outport
    allVars.outportIndex = zeros(initialSize,1); % index of specific outport. normally just 1 (will be more than 1 if the outport is a vector)
    allVars.optVarIdx = zeros(initialSize,1); % points to the true optimization variable. (will usually point to itself, otherwise, some "true" variable)
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
    block.inputs = cell(initialSize,1); % inputs of each block
    block.outportHandles = cell(initialSize,1); %all outport handles of block
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
            allVars,allVarsZero,block,blockZero,state,iZero,portH);
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
        
        %take into consideration if this block has been found before
        block.name{blockZero} = sourceBlock;
        block.handle(blockZero) = get_param(sourceBlock,'Handle');
        
        if strcmp(sourceType,'SubSystem')
            if isempty(refBlock)
                % if the current block is a subsystem and not from BLOM, 
                % want to look under the subsystem and get the appropriate
                % blocks there
                
                state = 'subsys';
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
                parentType = get_param(parentOfBlock,'BlockType');
                if strcmp(sourceType,'Inport') && strcmp(parentType,'SubSystem')
                    % in this case, there may actually be more inports
                    % connected
                    sourcePorts = get_param(parentOfBlock,'PortHandles');
                    sourceInports = [sourcePorts.Inport];
                    % storing the inport itself is not important. what we
                    % need to store is what is connected to the inport of
                    % the subsystem
                    outportHandles(iOut) = [];
                    iOut = iOut-1;
                    iZero = iZero-1;
                    if ~isempty(sourceInports)
                        [outportHandles,iZero,allVars,allVarsZero,block,blockZero] = ...
                            updateVars(sourceInports,outportHandles,iZero,...
                            allVars,allVarsZero,block,blockZero,2,iOut,sourcePorts);
                    else
                        iOut = iOut+1;
                        continue
                    end
                end
            else
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
    %check if any handles are -1 and remove them
    if any(outportHandles==-1)
        % FIX PRINT STATEMENT. NOT NECESSARILY TRUE
        fprintf('One or more of the blocks are missing a connection\n');
        stop = 1;
    end
    
    % remove empty and 0 entries in blocks
    for field={'names', 'P','K','inputs','outportHandles','dimensions',...
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
    end
    currentOutport = existingOutports(iOut);
    
    if strcmp(state,'subsys') %FIX: may want to look into special BLOM blocks here
        % if there's a subsystem, inports is actually an array of the
        % outports of subsystem
        parent = get_param(currentOutport,'Parent');
        index = inports==currentOutport;
        outportBlocks = find_system(parent,'regexp','on','BlockType','Outport');
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
        inports = [gotoPorts.Inport];
    end
    
    % found outports connected to inports provided
    for i = 1:length(inports);
        currentLine = get_param(inports(i),'Line');
        % this gives the all the outports connected to this line
        outportsFound = get_param(currentLine,'SrcPorthandle');
        outLength = length(outportsFound);
        % in case outportHandles is too short 
        if outLength > length(outportHandles)-iZero+1;
            if outlength > length(outportHandles)
                outportHandles = [outportHandles; zeros(outLength*2,1)];
            else
                outportHandles = [outportHandles; ...
                    zeros(length(outportHandles),1)];
            end
        end

        diff = setdiff(outportsFound,outportHandles);
        diffLength = length(diff);
        outportHandles(iZero:(diffLength+iZero-1)) = diff;
        iZero = iZero + diffLength;
    end
    
    switch state
        case 'bound'
            % bound case. for all the outports found here, include the
            % bounds for these
            allVarsState = 'bound';
            boundHandle = varargin{3};
            lowerBound = eval(get_param(boundHandle,'lb'));
            upperBound = eval(get_param(boundHandle,'ub'));
            
            % go through each of the outports connected through bound
            for currentOutport = outportsFound
                [block,blockZero,currentBlockIndex] =...
                    updateBlock(block,blockZero,currentOutport);
                [allVars,allVarsZero] = updateAllVars(allVars,allVarsZero,...
                    currentBlockIndex,currentOutport,allVarsState,lowerBound,upperBound);
            end
            
        case 'cost'
            
        otherwise 
            % not looking at bounds or costs
            allVarsState = 'normal';
            %update block structure
            [block,blockZero,currentBlockIndex] =...
                updateBlock(block,blockZero,currentOutport);
            [allVars,allVarsZero] = updateAllVars(allVars,allVarsZero,...
                    currentBlockIndex,currentOutport,allVarsState);
    end
        
    
end

%%
%======================================================================
%> @brief update allVars structure 
%>
%> More detailed description of the problem.
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

function [allVars,allVarsZero] = updateAllVars(allVars,allVarsZero,...
    currentBlockIndex,currentOutport,allVarsState,varargin)
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
                        'optVarIdx','cost'}
                        allVars.(field{1}) = [allVars.(field{1}); zeros(newLength*2,1)];
                end
                allVars.upperBound = [allVars.upperBound; inf*ones(newLength*2,1)];
                allVars.lowerBound = [allVars.lowerBound; -inf*ones(newLength*2,1)];
                allVars.time = [allVars.time; cell(newLength*2,1)];
            else % double the length of all fields in allVars
                for field={'block', 'outportNum','outportHandle','outportIndex',...
                        'optVarIdx','cost'}
                        allVars.(field{1}) = [allVars.(field{1}); zeros(oldLength,1)];
                end
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
                
            case 'normal'
                % do nothing
        end
        
        allVarsZero = allVarsZero+lengthOut;
    end
end

%%
%======================================================================
%> @brief given block structure and outport, populate the relevant fields
%> of block
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
%this function populates block using currentOutport information
    currentBlockHandle = get_param(currentOutport,'ParentHandle');
    referenceBlock = get_param(currentBlockHandle,'ReferenceBlock');
    % if the size of block equals the blockZero, need to double to
    % length of block
    if blockZero == length(block.handles)
        for field={'names', 'P','K','inputs','outportHandles','dimensions',...
                'sourceOutports'}
                block.(field{1}) = [block.(field{1}); cell(blockZero,1)];
        end
        block.handles = [block.handles; zeros(blockZero,1)];
    end

    if ~any(block.handles==currentBlockHandle)
        % no duplicate blocks, add this block
        currentBlockName = get_param(currentOutport,'Parent');
        sourcePorts = get_param(currentBlockName,'PortHandles');
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

        if isempty(block.outportHandles{blockZero})
            block.outportHandles{blockZero} = zeros(length(sourcePorts.Outport),1);
        end
        % FIX: POPULATE OUTPORT HANDLES

        % increase the index of block by one and populate currentBlockIndex
        currentBlockIndex = blockZero;
        blockZero = blockZero+1;
    else %case when the current outport's block is found but has not been added to block data
        currentBlockIndex = find(block.handles==currentBlockHandle);
        % FIX, find stuff for this case
    end
end

%%
%======================================================================
%> @brief Creates struct that says which outports are relevant at which
%> times
%>
%> More detailed description of the problem.
%>
%> @param outportHandles outportHandles found by searchSources
%>
%> @retval timeStruct structure with following fields. 1) outportHandles 2)
%> majorTimeStep 3) minorTimeStep
%======================================================================

function [timeStruct] = relevantTimes(outportHandles)
    timeStruct.outportHandles = outportHandles;
    % FIX: actually check where the time steps are relevant
    timeStruct.majorTimeStep = true(length(outportHandles),1);
    timeStruct.minorTimeStep = true(length(outportHandles),1);
end

%%
%======================================================================
%> @brief From outport list, create three different structures needed for
%> the optimization problem
%>
%> More detailed description of the problem.
%>
%> @param outportHandles outportHandles found by searchSources
%> @param name name of model
%>
%> @retval allVars structure with fields 1) index of block name 2) outport
%> # 3) outport handle 4) Index of outport 5) bounds 6) Cost 7) time 8)
%> list of indices for each optimization variable
%> @retval blocks structure with fields 1) name of block 2) handle of block
%> 3) P 4) K 5) inputs, allVars indices 6) outputs, optimvar indices 7)
%> all outport handles 8) dimensions 9) all source outport handles 
%======================================================================

function [allVars,polyStruct,blocks] = makeStruct(outportHandles,name)
    polyHandles = [find_system(name,'FindAll','on','ReferenceBlock',...
        'BLOM_Lib/Polyblock')];
    polyLength = length(polyHandles);
    
    % create structure to store P and K matricies
    polyStruct.block = cell(polyLength,1);
    polyStruct.P = cell(polyLength,1);
    polyStruct.K = cell(polyLength,1);
    polyIdx = 1;
    
    % structure for block names and handles
    blocks.names = cell(length(outportHandles),1);
    blocks.handles = zeros(length(outportHandles),1);
    blockZero = 1;
    
    % structure for optimization variables and handles
    allVars.block = zeros(length(outportHandles),1);
    allVars.index = zeros(length(outportHandles),1);
    allVarsZero = 1;
    for i = 1:length(outportHandles)
        currentBlockName = get_param(outportHandles(i),'Parent');
        currentBlockHandle = get_param(outportHandles(i),'ParentHandle');
        if ~any(blocks.handles==currentBlockHandle)
            % no duplicate blocks, add this block
            blocks.names{blockZero} = currentBlockName;
            blocks.handles(blockZero) = currentBlockHandle;
            blockZero = blockZero+1;
        end
        
        polyStruct.block{polyIdx} = blockZero-1;
        if any(polyHandles==currentBlockHandle)
            % store P&K matricies if the current block is a polyblock
            polyStruct.P{polyIdx} = eval(get_param(currentBlockHandle,'P'));
            polyStruct.K{polyIdx}= eval(get_param(currentBlockHandle,'K'));
        else
            % store P and K matricies for the other blocks
            [P,K] = BLOM_Convert2Polyblock(currentBlockHandle);
            polyStruct.P{polyIdx} = P;
            polyStruct.K{polyIdx}= K;
        end
        polyIdx = polyIdx +1;
        
        % get dimensions of each outport and for each parameter, store this
        % in a structure
        currentDim = get_param(outportHandles(i),'CompiledPortDimensions');
        lengthOut = currentDim(1)*currentDim(2);
        currentIndices = 1:lengthOut;
        if (allVarsZero+lengthOut) > length(allVars.block)
            % if the current length of optimvar is not enough, double the
            % size
            allVars.block = [allVars.block; zeros(length(allVars.block),1)];
            allVars.index = [allVars.index; zeros(length(allVars.index),1)];
        end
        % store the block and index
        allVars.block(allVarsZero:(allVarsZero+lengthOut-1)) = blockZero-1;
        allVars.index(allVarsZero:(allVarsZero+lengthOut-1)) = currentIndices;
        allVarsZero = allVarsZero+lengthOut;
        
    end
    
    blocks.names = blocks.names(1:blockZero-1);
    blocks.handles = blocks.handles(1:blockZero-1);
    allVars.block = allVars.block(1:allVarsZero-1);
    allVars.index = allVars.index(1:allVarsZero-1);
end

