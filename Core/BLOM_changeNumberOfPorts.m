%> @file BLOM_polyblockChangeInputs.m
%> @brief callback function for Polyblock that changes the number of inputs
%> to the block based on what the user specified
%>
%> More detailed description of the problem.
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
                block_name = [gcb '/In' num2str(i)];
                add_block('simulink/Sources/In1',block_name)
            end
        elseif inputNum < currentNumInputs
            for i = currentNumInputs:(inputNum+1)
                % delete unneccesary inports
                currentBlock = [gcb '/In' num2str(i)];
                delete_block(currentBlock);
            end
        else
            % do nothing. Number of inports equals number of desired inports.
            % This case should not be reached.
        end
    end
    
    if outputNum ~= -1
        currentNumOutputs = length(portInfo.Outport);
        if outputNum > currentNumOutputs
            for i = (currentNumOutputs+1):outputNum
                %create new number of needed outports
                block_name = [gcb 'Out' num2str(i)];
                add_block('simulink/Sources/Out1',block_name)
            end
        elseif outputNum < currentNumOutputs
            for i = currentNumOutputs:(outputNum+1)
                % delete unnessary outports
                currentBlock = [gcb '/Out' num2str(i)];
            end
        else
            % do nothing. Number of outports equals number of desired
            % outports. This should not be reached.
        end
    end
end