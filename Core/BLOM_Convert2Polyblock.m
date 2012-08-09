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
    inportPlaces.matrix
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
            % FIX: add support for 'Sum' block. the one block that has
            % default inputs on left and bottom and has a circle.
            inputsToAdd = get_param(blockHandle,'Inputs');
            subtract_indices = inputsToAdd=='-';
            if isempty(inportPlaces.matrix)
                % for n scalar inputs, P = eye(n), K = ones(1,n). for each ith
                % input that is subtracted, that index in K should equal -1
                % instead.  
                P = speye(totalInputs);
                K = ones(1,totalInputs); 

                K(subtract_indices) = -1;
            elseif length(inportPlaces.matrix) == 1 
                % one matrix/vector and the rest scalars
                P = eye(totalInputs);
                vectPlace = inportPlaces.matrix(1);
                vectLength = prod(inportDim{vectPlace});
                scalLength = length(inportPlaces.scalar);
                K = ones(vectLength,scalLength+vectLength);
                col = 1;
                for i = 1:length(inports)
                    if any(inportPlaces.scalar==i)
                        % current column is just the scalar input
                        if subtract_indices(i)
                            % current column is subtracted
                            K(:,col) = -1;
                        end
                        col = col+1;
                    else
                        %current column is the beginning of the vector
                        if subtract_indices(i)
                        % if the matrix/vector is being subtracted
                            K(:,(col:(vectLength+col-1))) = -1*eye(vectLength);
                        else
                        % if it's being added
                            K(:,(col:(vectLength+col-1))) = eye(vectLength);
                        end
                        col = col+vectLength;
                    end
                end
            elseif isempty(inportPlaces.scalar)
                % all matricies/vectors. assumes that they are all the same
                % size. add element wise
                vectPlace = inportPlaces.matrix(1);
                vectLength = prod(inportDim{vectPlace});
                numVect = length(inports);
                P = speye(vectLength*numVect);
                K = zeros(vectLength,vectLength*numVect);
                j = 1;
                for i = 1:vectLength:(numVect*vectLength)
                    if subtract_indices(j)
                        K(:,(i:(i+vectLength-1))) = -1*eye(vectLength);
                    else
                        K(:,(i:(i+vectLength-1))) = eye(vectLength);
                    end
                    j = j+1;
                end
            elseif ~isempty(inportPlaces.scalar) && ~isempty(inportPlaces.matrix)
                % more than 1 matrix/vector of the same size and one or
                % more scalars
                vectPlace = inportPlaces.matrix(1);
                vectLength = prod(inportDim{vectPlace});
                numVect = length(inportPlaces.matrix);
                P = speye(vectLength*numVect+length(inportPlaces.scalar));
                K = ones(vectLength,vectLength*numVect+length(inportPlaces.scalar));
                col = 1;
                for i = 1:length(inports)
                    if any(inportPlaces.scalar==i)
                        % current column is just the scalar input
                        if subtract_indices(i)
                            % current column is subtracted
                            K(:,col) = -1;
                        end
                        col = col+1;
                    else
                        %current column is the beginning of the vector
                        if subtract_indices(i)
                        % if the matrix/vector is being subtracted
                            K(:,(col:(vectLength+col-1))) = -1*eye(vectLength);
                        else
                        % if it's being added
                            K(:,(col:(vectLength+col-1))) = eye(vectLength);
                        end
                        col = col+vectLength;
                    end
                end
            else
                % most likely will not reach this case
                P = [];
                K = [];
            end
        %% absolute value (only for costs)    
        case 'Abs'
            
            
        %% multiply/divide. Divide is the same block   
        case 'Product' 
            % for n scalar inputs, P = ones(1,n), K = [1]. for each ith
            % input to be subtracted, that index in P should be equal -1
            
            % only_mult returns 1 if it is just multiplication. Returns 0 if
            % there is division too
            inputs = get_param(blockHandle,'Inputs');
            digitInput = isstrprop(inputs,'digit');
            division = any(digitInput==0);
            if division
                divide_indices = inputs=='/';
            end
            
            mult_type = get_param(blockHandle,'Multiplication');
            switch mult_type
                case 'Element-wise(.*)'
                    % element by element multiplication
                    if isempty(inportPlaces.matrix)
                        % all scalars
                        P = ones(1,totalInputs);
                        K = 1;
                        if division
                            P(divide_indices) = -1;
                        end
                    elseif isempty(inportPlaces.scalar)
                        % all vectors/matrices of same size. multiply
                        % element wise
                        vectPlace = inportPlaces.matrix(1);
                        vectLength = prod(inportDim{vectPlace});
                        numVect = length(inports);
                        P = zeros(vectLength,vectLength*numVect);
                        K = speye(vectLength);
                        divide_index = 1;
                        for i = 1:vectLength:(numVect*vectLength)
                            if division
                                if divide_indices(divide_index)
                                    P(:,(i:(i+vectLength-1))) = -1*eye(vectLength);
                                else
                                    P(:,(i:(i+vectLength-1))) = eye(vectLength);
                                end
                            else
                                P(:,(i:(i+vectLength-1))) = eye(vectLength);
                            end
                            divide_index = divide_index+1;
                        end
                    elseif ~isempty(inportPlaces.scalar) && ~isempty(inportPlaces.matrix)
                        % one more more vectors/matrices and one more more
                        % scalars. all the vectors/matrices are multiplied
                        % element wise and all elements of the
                        % vectors/matrices are multiplied by the scalars
                        vectPlace = inportPlaces.matrix(1);
                        vectLength = prod(inportDim{vectPlace});
                        K = speye(vectLength);
                        P = ones(vectLength,vectLength*length(inportPlaces.matrix)...
                            +length(inportPlaces.scalar));
                        col = 1;
                        divide_index = 1;
                        for i = 1:length(inports)
                            if any(inportPlaces.matrix==i)
                                % current column is the vector
                                if division
                                    if divide_indices(divide_index)
                                        P(:,(col:(vectLength+col-1))) = -1*eye(vectLength);
                                    else
                                        P(:,(col:(vectLength+col-1))) = eye(vectLength);
                                    end
                                    col = col+vectLength;
                                    divide_index = divide_index+1;
                                else
                                    P(:,(col:(vectLength+col-1))) = eye(vectLength);
                                    col = col+vectLength;
                                end
                            else
                                if division
                                    if divide_indices(divide_index)
                                        P(:,col) = -1;
                                    end
                                    divide_index=divide_index+1;
                                end
                                col = col+1;
                            end
                        end
                        P = sparse(P);
                    end
                case 'Matrix(*)'
                    P = [];
                    K = [];
                    fprintf('Matrix multiplication currently not supported by BLOM\n');
            end
        %% constant    
        case 'Constant'
            if isempty(inportPlaces.matrix)
                P = 1;
                % FIX: find right value for K
                K = [];
            else
                vectPlace = inportPlaces.matrix(1);
                vectLength = prod(inportDim{vectPlace});
                P = speye(vectLength);
                % FIX: find right value for K
                K = [];
            end
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
    % gives indices of scalars and vector inports from dimensions
    % totalInputs gives the total number of elements (each element of a
    % matrix or vector adds to this number)
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

function [P_real,K_real] = friendlyToReal(P,K,outportDim)
    % converts the user-friendly P and K matrices to the proper P and K
    % matrices
    
    total_outputs = 0;
    % first figure out the number of outputs
    for i = 1:length(outportDim)
        total_outputs = total_outputs + prod(outportDim{i});
    end
    
    [M_P,N_P] = size(P);
    [M_K,N_K] = size(K);
    
    P_real = sparse(M_P+total_outputs,N_P+total_outputs);
    K_real = sparse(M_K,N_P+total_outputs);
    
    P_real(1:M_P,1:N_P) = P;
    K_real(1:M_K,1:N_K) = K;
    
    P_real( (M_P+1):end , (N_P+1):end ) = speye(total_outputs);
    K_real( :, (N_P+1):end) = -1*speye(total_outputs);
end