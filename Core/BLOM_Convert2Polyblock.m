%> @file BLOM_convert2Polyblock.m
%> @brief function that converts regular simulink blocks to a polyblock
%> like code. it will output the P&K matricies for a given Simulink block
%>
%> More detailed description of the problem.
%>
%> @param blockHandle handle of the Simulink block to be converted
%> @param varargin the external and input from simulink handles. the
%> sources of these blocks are not relevant to the optimization problem
%>
%> @retval P P matrix for f
%> @retval K K matrix
%======================================================================

function [P,K] = BLOM_Convert2Polyblock(blockHandle)
    % figure out block type
    blockType = get_param(blockHandle,'BlockType');
    portHandles = get_param(blockHandle,'PortHandles');
    inports = portHandles.Inport;
    outports = portHandles.Outport;

    % FIX: not sure if getting all the dimensions in a cell array is faster
    % than using just an array. Check to see which is faster
    inportDim = get_param(inports,'CompiledPortDimensions');
    outportDim = get_param(outports,'CompiledPortDimensions');
    
    % give port indices of vectors and scalars for inports
    [inportPlaces,totalInputs] = scalarVectIndex(inportDim);
    
%     inportDim = array(length(inports),2);
%     outportDim = array(length(outports),2);
%     
%     for i = 1:length(inports)
%         inportDim(i,:) = get_param(inports(i),'CompiledPortDimensions');
%     end
%     
%     for i = 1:length(outports)
%         outportDim(i,:) = get_param(outports(i),'CompiledPortDimensions');
%     end
    
    %% switch for different cases. hopefully able to compare numbers
    %% instead of strings
    switch blockType
        case {'Sum','Add'} % add/subtract
            if isempty(inportPlaces.matrix)
                % for n scalar inputs, P = eye(n), K = ones(1,n). for each ith
                % input that is subtracted, that index in K should equal -1
                % instead.  
                portSize = length(totalInputs);
                P = speye(portSize(1));
                K = ones(1,portSize(1)); 

                inputsToAdd = get_param(blockHandle,'Inputs');
                subtract_indicies = inputsToAdd=='-';
                K(subtract_indicies) = -1;
            elseif length(inportPlaces.matrix) == 1 
                % one matrix/vector and the rest scalars
                P = eye(totalInputs);
                vectPlace = inportPlaces.matrix(1);
                vectLength = prod(inportDim{vectPlace});
                scalLength = length(inportPlaces.scalar);
                K = ones(vectLength,scalLength+vectLength);
                K(:,(vectPlace:(vectLength+vectPlace-1))) = eye(vectLength);
            end
        %% absolute value (only for costs)    
        case 'Abs'
            
            
        %% multiply/divide     
        case {'Product','Divide'}
            % for n scalar inputs, P = ones(1,n), K = [1]. for each ith
            % input to be subtracted, that index in P should be equal -1
            
        %% constant    
        case 'Constant'
            
        %% gain
        case 'Gain'
            
        %% bias
        case 'Bias'
            
        %% trigonometric function (currently not functional in polyblock)
        case 'Trigonometric Function'
            
        otherwise 
            P = [];
            K = [];
            fprintf('This block is currently not supported by BLOM\n');

    end


end

function [inportPlaces,totalInputs] = scalarVectIndex(dimCell)
    % gives indicies of scalars and vector inports from dimensions
    inportPlaces.scalar = zeros(length(dimCell),1);
    inportPlaces.matrix = zeros(length(dimCell),1);
    scalarZero = 1;
    matrixZero = 1;
    totalInputs = 0;
    for i = 1:length(dimCell);
        if prod(dimCell{i}) == 1
            inportPlaces.scalar(scalarZero) = i;
            scalarZero = scalarZero+1;
            totalInputs = totalInputs+1;
        else
            inportPlaces.matrix(matrixZero) = i;
            matrixZero = matrixZero+1;
            totalInputs = totalInputs + prod(dimCell{i});
        end
    end
    
    if scalarZero ~= 1
        inportPlaces.scalar = inportPlaces.scalar(1:(scalarZero-1));
    else
        inportPlaces.scalar = [];
    end
    
    if matrixZero ~= 1
        inportPlaces.matrix = inportPlaces.matrix(1:(matrixZero-1));
    else
        inportPlaces.matrix = [];
    end
    
end