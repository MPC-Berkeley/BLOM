%> @file BLOM_polyblockChangeInputs.m
%> @brief callback function for Polyblock that changes the number of inputs
%> to the block based on what the user specified
%>
%> In it's current state, this program assumes that all the inputs and
%> outputs of the BLOM Library blocks have a standard naming convention,
%> e.g. In1, In2,... and Out1, Out2...
%>
%> @param block usually gcb
%> @param inputNum number of desired inputs. -1 if no change
%> @param outputNum number of desired outputs. -1 if no change
%>
%======================================================================

function BLOM_changeNumberOfPorts(block,inputNum,outputNum)
    portInfo = get_param(block,'PortHandles');
    
    if inputNum ~= -1
        % find number of inports
        currentNumInputs = length(portInfo.Inport);
        if inputNum > currentNumInputs
            for i = (currentNumInputs+1):inputNum
                % create new number of needed inports
                block_name = [block '/In' num2str(i)];
                add_block('simulink/Sources/In1',block_name)
            end
        elseif inputNum < currentNumInputs
            for i = currentNumInputs:-1:(inputNum+1)
                % delete unneccesary inports
                currentBlock = [block '/In' num2str(i)];
                delete_block(currentBlock);
            end
        else
            % do nothing. Number of inports equals number of desired
            % inports.
        end
    end
    
    if outputNum ~= -1
        currentNumOutputs = length(portInfo.Outport);
        if outputNum > currentNumOutputs
            for i = (currentNumOutputs+1):outputNum
                %create new number of needed outports
                block_name = [block '/Out' num2str(i)];
                add_block('simulink/Sinks/Out1',block_name)
            end
        elseif outputNum < currentNumOutputs
            for i = currentNumOutputs:-1:(outputNum+1)
                % delete unnessary outports
                currentBlock = [block '/Out' num2str(i)];
                delete_block(currentBlock);
            end
        else
            % do nothing. Number of outports equals number of desired
            % outports.
        end
    end
end