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
    
    % figure out total number of outputs
    total_outputs = 0;
    for i = 1:length(outportDim)
        class_outportDim = outportDim;
        if iscell(class_outportDim)
            total_outputs = total_outputs + prod(outportDim{i});
        else
            total_outputs = prod(outportDim);
        end
    end
    
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

            % go from user friendly form of P and K to proper P and K
            % matricies
            P = blkdiag(P,speye(total_outputs));
            K = horzcat(K,-1*speye(total_outputs));
            
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
                        % one or more vectors/matrices and one or more
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
            
            % go from user friendly form of P and K to proper P and K
            % matricies
            P = blkdiag(P,speye(total_outputs));
            K = horzcat(K,-1*speye(total_outputs));
        %% constant    
        case 'Constant'
            P = vertcat(speye(total_outputs),zeros(1,total_outputs));
            K = horzcat(-1*speye(total_outputs),zeros(total_outputs,1));
            constVal = get_param(blockHandle,'Value');
            constVal = evalin('base',constVal);
            K(:,end) = constVal;
        %% gain
        case 'Gain'
            % FIX, currently no support for matrix multiplication
            mult_type = get_param(blockHandle,'Multiplication');
            switch mult_type
                case 'Element-wise(K.*u)'
                    P = speye(totalInputs*2);
                    gain = evalin('base',get_param(blockHandle,'Gain'));
                    K = horzcat(gain*speye(totalInputs),-1*speye(totalInputs));
                otherwise
                    P = [];
                    K = [];
                    fprintf('Matrix multiplication currently not supported by BLOM\n');
            end
        %% bias
        case 'Bias'
            P = speye(totalInputs*2+1);
            P(end) = 0;
            
            bias = evalin('base',get_param(blockHandle,'Bias'));
            
            K = horzcat(speye(totalInputs),-1*eye(totalInputs),...
                bias*ones(totalInputs,1));
        %% trigonometric function 
        
        case 'Trigonometry'
            % expresses all trig functions using a unique value e20 (e.g. sin
            % is 3e20)   
            mathFunction=get_param(blockHandle,'Operator');
            P=blkdiag(speye(totalInputs)*BLOM_FunctionCode(mathFunction), speye(totalInputs));
            K=horzcat(speye(totalInputs),-1*speye(totalInputs));
        %% polynomial    
        case 'Polyval'
            polyInportDim=get_param(blockHandle,'CompiledPortDimensions');
            % length of input vector
            polyInputDimNum=prod(polyInportDim.Inport);
            % get coefficients
            coeffs=eval(get_param(blockHandle,'Coefs'));
            % number of coefficients
            coeffsLength=length(coeffs); 
            % vector of powers of polynomial
            polyPowers=linspace(coeffsLength-1,0,coeffsLength)'; 
            % vector of power of y (always zeros then a 1)
            polyPowersy=[zeros(coeffsLength-1,1);1]; 
            % left half of P matrix
            pLeft=zeros(coeffsLength*polyInputDimNum,polyInputDimNum); 
            % right half of P matrix
            pRight=pLeft; 
            K=zeros(polyInputDimNum,polyInputDimNum*coeffsLength);
            
            % create left and right halves of P matrix, as well as K matrix
            for polyIndex=1:polyInputDimNum 
                pLeft(((polyIndex-1)*coeffsLength+1):polyIndex*coeffsLength,polyIndex)=polyPowers;
                pRight(((polyIndex-1)*coeffsLength+1):polyIndex*coeffsLength,polyIndex)=polyPowersy;
                K(polyIndex,(polyIndex-1)*coeffsLength+1:polyIndex*coeffsLength)=[coeffs(1:end-1) -1];
            end
            % combine left and right halves
            P=[pLeft pRight]; 
            % insert elements corresponding to constants at the end
            P=sparse([P;zeros(1,polyInputDimNum*2)]); 
            % insert elements corresponding to constants at the ends
            K=sparse([K coeffs(end)*ones(polyInputDimNum,1)]); 
            
        %% unary minus
        case 'UnaryMinus'
                    P = speye(totalInputs*2);
                    K = horzcat(-speye(totalInputs),-1*speye(totalInputs));
                    
        %% Math Function
        case 'Math'
            % Currently expresses exp, log, log10, conj, rem,
            % mod using a unique value e20 for each special function. The
            % other special functions will be expressed using P and K
            % matrices.
            mathFunction=get_param(blockHandle,'Operator');
            P=blkdiag(speye(totalInputs)*BLOM_FunctionCode(mathFunction), speye(totalInputs));
            K=horzcat(speye(totalInputs),-1*speye(totalInputs));
            
        %% SubSystem.
        case 'SubSystem'
            % currently just handles the bound blom block. Other blocks
            % with subsystem block types return an empty P and K
            
            % find the type of the block (Bound, DiscreteCost, etc)
            referenceBlock=get_param(blockHandle,'ReferenceBlock'); 
            blomBlockType=referenceBlock(10:end);
            
            % bound block uses inequality constraints, no bound exists if
            % Inf or -Inf is present             
            if strcmp(blomBlockType,'Bound')==1
                % extract upper and lower bounds
                upperbound=eval(get_param(blockHandle,'ub'));
                lowerbound=eval(get_param(blockHandle,'lb'));
                P=sparse([ones(2*totalInputs,1) zeros(2*totalInputs,1)]);
                K=sparse([-ones(totalInputs,1) lowerbound*ones(totalInputs,1);ones(totalInputs,1) -upperbound*ones(totalInputs,1)]);
                if upperbound==Inf
                    P(totalInputs+1:end,:)=[];
                    K(totalInputs+1:end,:)=[];
                end
                if lowerbound==-Inf
                    P(1:totalInputs,:)=[];
                    K(1:totalInputs,:)=[];
                end
            else
                P=[];
                K=[];
            end
            
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
    if iscell(dimCell)
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
    else
        if prod(dimCell) == 1
            inportPlaces.scalar(scalarZero) = 1;
            scalarZero = scalarZero+1;
            totalInputs = totalInputs+1;
        else
            inportPlaces.matrix(matrixZero) = 1;
            matrixZero = matrixZero+1;
            totalInputs = totalInputs + prod(dimCell);
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