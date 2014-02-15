%> @file BLOM_Subsys2Poly.m
%> @brief function that allows you to convert a simple subsystem into a
%> polyblock by providing you with P,K of internal contents.  To call on
%> group of blocks, highlight and make into subsystem then call this
%> function on resulting subsystem
%>
%> @param subsysHandle handle of subsystem to be converted
%>
%> @retval P matrix to input into polyblock
%> @retval K matrix to input into polyblock
%========================================================================
function [P K] = BLOM_Subsys2Poly(subsysHandle)
    warning('Currently does not support vectors');
    
    outportBlocks = find_system(subsysHandle,'SearchDepth',1,'regexp','on','BlockType','Outport');
    portCount = get_param(subsysHandle, 'Ports');
    numInports = portCount(1);

    P = [];
    K = [];

    for idx = 1:length(outportBlocks)
        outportBlock = outportBlocks(idx);
        ports = get_param(outportBlock, 'PortHandles');
        inport = ports.Inport;
        [newP newK] = getPK(inport, numInports);
        P = [P; newP];
        K = blkdiag(K, newK);
    end
    
    [P K] = cleanupPK(P, K);

    %=====================================================%
    % uncomment below to auomatically replace subsystem   %
    % WARNING: make sure to save before trying to replace %
    % it will delete your previous blocks and you might   %
    % not be able to undo.                                %
    %=====================================================%
    
%     [Pm, Pn, Ps] = find(P);
%     if size(Pm,1) == 1
%         Pstr = ['sparse([' num2str(Pm) '],[' num2str(Pn) '],[' num2str(Ps) '])'];
%     elseif size(Pm,2)==1
%         Pstr = ['sparse([' num2str(Pm') '],[' num2str(Pn') '],[' num2str(Ps') '])'];
%     else
%         error('P dimensions are weird... should never get here');
%     end
%     [Km, Kn, Ks] = find(K);
%     if size(Km,1)==1
%         Kstr = ['sparse([' num2str(Km) '],[' num2str(Kn) '],[' num2str(Ks) '])'];
%     elseif size(Km,2)==1
%         Kstr = ['sparse([' num2str(Km') '],[' num2str(Kn') '],[' num2str(Ks') '])'];
%     else
%         error('K dimensions are weird... should never get here');
%     end
%     
%     name = get_param(subsysHandle, 'Name');
%     parent = get_param(subsysHandle, 'Parent');
%     position = get_param(subsysHandle, 'Position');
%     orientation = get_param(subsysHandle, 'Orientation');
%     block = [parent '/' name];
%     
%     delete_block(subsysHandle);
%     add_block('BLOM_Lib/Polyblock', block);
%     
%     set_param(block, 'P', Pstr);
%     set_param(block, 'K', Kstr);
%     set_param(block, 'inputScalar', 'on');
%     set_param(block, 'outputScalar', 'on');
%     set_param(block, 'Orientation', orientation);
%     set_param(block, 'Position', position);    
    
end

%%
%======================================================================
%> @brief gives PK matrix starting from given inport feeding into a block
%>
%> @param portHandle to start from
%> @param numInputs number of inports to subsystem (also number of input
%> variables)
%>
%> @retval P matrix
%> @retval K matrix
%=======================================================================
function [P K] = getPK(portHandle, numInputs)
    line = get_param(portHandle, 'Line');
    parent = get_param(line, 'SrcBlockHandle');
    blockType = get_param(parent, 'BlockType');
    blockName = get_param(parent, 'Name');
    
    switch blockType
        case 'Inport'
            inputNum = eval(get_param(parent, 'Port'));
            P = sparse(1, inputNum, 1, 1, numInputs);
            K = sparse(1,1,1);
        
        case 'Constant'
            P = spalloc(1,numInputs,numInputs);
            K = sparse(1,1,eval(get_param(parent, 'Value')));
        
        case {'Sum', 'Add'}
            inputSigns = get_param(parent, 'Inputs');
            inputSigns(inputSigns == '|') = [];
            ports = get_param(parent, 'PortHandles');
            inports = ports.Inport;
            P = [];
            K = [];
            for input = 1:length(inports)
                [newP newK] = getPK(inports(input), numInputs);
                P = [P;newP];
                if inputSigns(input) == '+'
                    K = [K newK];
                elseif inputSigns(input) == '-'
                    K = [K -newK];
                else
                    error(['Invalid input sign ''' inputSigns(input) ''' in block ' blockName]);
                end                    
            end            
            
        case 'UnaryMinus'
            ports = get_param(parent, 'PortHandles');
            inport = ports.Inport;
            [P newK] = getPK(inport, numInputs);
            K = -newK;
            
        case 'Gain'
            gain = eval(get_param(parent, 'Gain'));
            ports = get_param(parent, 'PortHandles');
            inport = ports.Inport;
            [P newK] = getPK(inport, numInputs);
            K = gain * newK;
            
        case 'Bias'
            bias = eval(get_param(parent, 'Bias'));
            ports = get_param(parent, 'PortHandles');
            inport = ports.Inport;
            [newP newK] = getPK(inport, numInputs);
            P = [newP; spalloc(1,numInputs,numInputs)];
            K = [newK bias];
            
        case 'From'
            tag = get_param(parent, 'GoToTag');
            system = get_param(parent, 'Parent');
            gotoBlock = find_system(system, 'BlockType', 'Goto', 'GotoTag', tag);
            gotoPorts = get_param(gotoBlock{1},'PortHandles');
            inport = [gotoPorts.Inport];
            [P K] = getPK(inport, numInputs);
            
        case 'Product'           
            inputs = get_param(parent, 'Inputs');
            division = any(inputs == '/');
            if division
                error(['Invalid block: ' blockName '. Currently do not support division'])
            end
            ports = get_param(parent, 'PortHandles');
            inports = ports.Inport;
            P = spalloc(1,numInputs,numInputs);
            K = 1;
            
            for input = 1:length(inports)
                [newP newK] = getPK(inports(input), numInputs);
                P = repmat(P,size(newP,1),1) + kron(newP, ones(size(P,1),1));
                Kmat = K' * newK;
                K = Kmat(1:end);            
            end 
            
        otherwise
            error(['Cannot create Polyblock.  ' blockName ' is of unsupported type: ' blockType]);
    end
end

%%
%======================================================================
%> @brief removes redundant rows of P and combines corresponding columns of
%> K and remove all empty columns of K and corresponding rows of P
%>
%> @param P matrix to cleanup
%> @param K matrix to cleanup
%>
%> @retval newP matrix reduced
%> @retval newK matrix reduced
%=======================================================================
function [newP newK] = cleanupPK(P, K)
    [newP,~,J] = unique(P, 'rows', 'last');
    newK = zeros(size(K,1), max(J));
    for col = 1:size(K,2)
       newK(:,J(col)) =  newK(:,J(col)) + K(:,col);
    end

    removeTerms = find(~any(newK));
    newK(:,removeTerms) = [];
    newP(removeTerms,:) = [];
end