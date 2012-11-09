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

    [boundHandles,costHandles,inputAndExternalHandles] = findBlocks(name);
    [outportHandles,boundStruct,block,optimVar,stop] = ...
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
    fprintf('The Number of outports is %.0f\n',length(outportHandles))
    for i = 1:length(outportHandles);
        parent = get_param(boundStruct.outportHandles(i),'Parent')
        portType = get_param(outportHandles(i),'PortType');
    end
    
    fprintf('\n\n\n Here we have the stored data from block\n\n\n');
    
    for i = 1:length(block.handles);
        someBlock = block.names{i}
    end
    
    % just a placeholder for ModelSpec so that MATLAB does not complain
    ModelSpec = 1;
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

function [outportHandles,boundStruct,block,optimVar,stop] = ...
    searchSources(boundHandles,costHandles,varargin)
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
    
    optimZero = 1; % index of first optimVar zero
    % optimVar stores information for every optimization variable
    optimVar.block = zeros(initialSize,1); % the index of the block
    optimVar.outportNum = zeros(initialSize,1); % outport number
    optimVar.outportHandle = zeros(initialSize,1); % handle of specific outport
    optimVar.outportIndex = zeros(initialSize,1); % index of specific outport. normally just 1
    optimVar.bounds = cell(initialSize,1); % store the bounds of each variable. will be an array stored in each cell element
    optimVar.cost = cell(initialSize,1);
    optimVar.time = cell(initialSize,1);
    optimVar.sameOpt = cell(initialSize,1); % a list of indices for duplicate variables
    
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
    % ports from there. fill in optimVar and block structures as necessary
    
    % iZero is index of first zero
    iZero = 1;
    outportHandles = zeros(initialSize,1);    
    % get all outports connected to the bounds 
    for i = 1:length(boundHandles)
        portH = get_param(boundHandles(i),'PortHandles');
        % costs and bounds should only have one inport and line
        currentInport = portH.Inport;
        [outportHandles,iZero,optimVar,optimZero,block,blockZero] = ...
            updateVars(currentInport,outportHandles,iZero,...
            optimVar,optimZero,block,blockZero,1,iZero,portH);
    end
    
    % create structure for bounds. fill entries up to iBounds with true.
    % These are all the outportHandles found so far and must have a bound
    iBounds = iZero-1;
    boundStruct.bound = true(iBounds,1);
    
    % get all outports connected to costs
    for i = 1:length(costHandles)
        portH = get_param(costHandles(i),'PortHandles');
        % costs and bounds should only have one inport and line
        currentInport = portH.Inport;
        [outportHandles,iZero,optimVar,optimZero,block,blockZero] = ...
            updateVars(currentInport,outportHandles,iZero,...
            optimVar,optimZero,block,blockZero,1,iZero,portH);
    end
    
    % array of outport indices to remove at the end of search. for now, we
    % want to remove the outports of From blocks. Do not want to remove it
    % during the search to prevent finding it over again each time
    removeOutport = zeros(10,1);
    removeOutportZeroIndex = 1;
    
    
    % iOut is current handle we are looking at
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
        
        if strcmp(sourceType,'SubSystem') && isempty(refBlock)
            % if the current block is a subsystem and not from BLOM, 
            % want to look under the subsystem and get the appropriate
            % blocks there
            sourceOutports = [sourcePorts.Outport];
            [outportHandles,iZero,optimVar,optimZero,block,blockZero] = ...
                updateVars(sourceOutports,outportHandles,iZero,...
                optimVar,optimZero,block,blockZero,1,iOut,sourcePorts);
            
            % want to remove this outport later since all is does is route
            % a signal
            removeOutport(removeOutportZeroIndex) = iOut;
            removeOutportZeroIndex = removeOutportZeroIndex + 1;
            
            if removeOutportZeroIndex==length(removeOutport)
                removeOutport = [removeOutport; zeros(length(removeOutport),1)];
            end
        elseif strcmp(sourceType,'From')
            % the current block is a from block, find goto block and see
            % what is connected to the goto block
            tag = get_param(sourceBlock,'GotoTag');
            
            % want to remove this outport later since all is does is route
            % a signal
            removeOutport(removeOutportZeroIndex) = iOut;
            removeOutportZeroIndex = removeOutportZeroIndex + 1;
            
            if removeOutportZeroIndex==length(removeOutport)
                removeOutport = [removeOutport; zeros(length(removeOutport),1)];
            end
            
            % there should only be one goto block associated with this from
            % block
            gotoBlock = find_system(name,'BlockType','Goto','GotoTag',tag);
            gotoPorts = get_param(gotoBlock{1},'PortHandles');
            sourceInports = [gotoPorts.Inport];
            [outportHandles,iZero,optimVar,optimZero,block,blockZero] = ...
                updateVars(sourceInports,outportHandles,iZero,...
                optimVar,optimZero,block,blockZero,1,iOut,sourcePorts);
            
        elseif length(sourceInports) > 1
            % currently do nothing special if there is more than one input
            % except look for what's connected. 
            % FIX: check to see which inports affect which outports. Not
            % all inports may be relevant
            [outportHandles,iZero,optimVar,optimZero,block,blockZero] = ...
                updateVars(sourceInports,outportHandles,iZero,...
                optimVar,optimZero,block,blockZero,2,iOut,sourcePorts);
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
                        [outportHandles,iZero,optimVar,optimZero,block,blockZero] = ...
                            updateVars(sourceInports,outportHandles,iZero,...
                            optimVar,optimZero,block,blockZero,2,iOut,sourcePorts);
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
            [outportHandles,iZero,optimVar,optimZero,block,blockZero] = ...
                updateVars(sourceInports,outportHandles,iZero,...
                optimVar,optimZero,block,blockZero,2,iOut,sourcePorts);
        end
        % FIX: should check to see if BLOM supports the blocks that is
        % found right here
        iOut = iOut+1;
    end
    
    % remove unwanted outports (e.g. from block outport)
    if any(removeOutport==0)
        removeOutport = removeOutport(1:(removeOutportZeroIndex-1));
    end
    outportHandles(removeOutport) = [];
    iZero = iZero-removeOutportZeroIndex+1;
    
    if any(outportHandles==0)
        outportHandles = outportHandles(1:iZero-1);
    end
    %check if any handles are -1 and remove them
    if any(outportHandles==-1)
        % FIX PRINT STATEMENT. NOT NECESSARILY TRUE
        fprintf('One or more of the blocks are missing a connection\n');
        stop = 1;
    end
    
    % FIX: need to find some way to remove all -1 handles. using setdiff
    % with [-1] reorders all the outport handles and puts it in ascending
    % order
    % outportHandles = setdiff(outportHandles,[-1]);
    
    % update boundStruct with all the outport handles and boolean values
    % for true and false
    boundStruct.bound = [boundStruct.bound; false(length(outportHandles)-iBounds,1)];
    boundStruct.outportHandles = outportHandles;
end

%======================================================================
%> @brief From given inports, see which outports are relevant. Also, given
%> the current outport, populate the block and optimVar structures
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

function [outportHandles,iZero,optimVar,optimZero,block,blockZero] =...
    updateVars(inports,existingOutports,iZero,optimVar,optimZero,...
    block,blockZero,state,varargin)

    outportHandles = existingOutports;
    % if the current state is a subsystem, do special case for subsystem
    if state == 3
        % the current block is a subsystem
        subsys = true;
    else
        subsys = false;
    end
    
    if ~isempty(varargin)
        iOut = varargin{1};
        sourcePorts = varargin{2};
    end
    
    if state ~= 1 % not looking at bounds or costs
        currentOutport = existingOutports(iOut);
        currentBlockHandle = get_param(currentOutport,'ParentHandle');
        referenceBlock = get_param(currentBlockHandle,'ReferenceBlock');

        % if the size of block equals the blockZero, need to double to
        % length of block
        if blockZero == length(block.handles)
            for field={'name', 'P','K','inputs','outportHandles','dimensions',...
                    'sourceOutports'}
                    block.(field{1}) = [strarr.(field{1}); cell(blockZero,1)];
            end
            block.handles = [block.handles; zeros(blockZero,1)];
        end

        if ~any(block.handles==currentBlockHandle)
            % no duplicate blocks, add this block
            currentBlockName = get_param(currentOutport,'Parent');
            block.names{blockZero} = currentBlockName
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

            % increase the index of block by one
            blockZero = blockZero+1;
        end
    end
    
    
    switch state
        case 3 % current block is a subsystem
            currentOutport = existingOutports(iOut);
            % if there's a subsystem, inports is actually an array of the
            % outports
            parent = get_param(currentOutport,'Parent');
            index = inports==currentOutport;
            outportBlocks = find_system(parent,'regexp','on','BlockType','Outport');
            handle = get_param(outportBlocks{index},'Handle');
            portH = get_param(handle,'PortHandles');
            currentInport = [portH.Inport];
            [outportHandles,iZero] = getOutports(currentInport,outportHandles,...
                iZero);
        otherwise
            for i = 1:length(inports);
                currentLine = get_param(inports(i),'Line');
                % this gives the all the outports connected to this line
                currentOutports = get_param(currentLine,'SrcPorthandle');
                outLength = length(currentOutports);
                % in case outportHandles is too short 
                if outLength > length(outportHandles)-iZero+1;
                    if outlength > length(outportHandles)
                        outportHandles = [outportHandles; zeros(outLength*2,1)];
                    else
                        outportHandles = [outportHandles; ...
                            zeros(length(outportHandles),1)];
                    end
                end

                diff = setdiff(currentOutports,outportHandles);
                diffLength = length(diff);
                outportHandles(iZero:(diffLength+iZero-1)) = diff;
                iZero = iZero + diffLength;
            end
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
%> @retval optimVar structure with fields 1) index of block name 2) outport
%> # 3) outport handle 4) Index of outport 5) bounds 6) Cost 7) time 8)
%> list of indices for each optimization variable
%> @retval blocks structure with fields 1) name of block 2) handle of block
%> 3) P 4) K 5) inputs, optimVar indices 6) outputs, optimvar indices 7)
%> all outport handles 8) dimensions 9) all source outport handles 
%======================================================================

function [optimVar,polyStruct,blocks] = makeStruct(outportHandles,name)
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
    optimVar.block = zeros(length(outportHandles),1);
    optimVar.index = zeros(length(outportHandles),1);
    optimZero = 1;
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
        if (optimZero+lengthOut) > length(optimVar.block)
            % if the current length of optimvar is not enough, double the
            % size
            optimVar.block = [optimVar.block; zeros(length(optimVar.block),1)];
            optimVar.index = [optimVar.index; zeros(length(optimVar.index),1)];
        end
        % store the block and index
        optimVar.block(optimZero:(optimZero+lengthOut-1)) = blockZero-1;
        optimVar.index(optimZero:(optimZero+lengthOut-1)) = currentIndices;
        optimZero = optimZero+lengthOut;
        
    end
    
    blocks.names = blocks.names(1:blockZero-1);
    blocks.handles = blocks.handles(1:blockZero-1);
    optimVar.block = optimVar.block(1:optimZero-1);
    optimVar.index = optimVar.index(1:optimZero-1);
end

