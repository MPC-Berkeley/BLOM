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
        muxSource = [block '/Mux'];
        if inputNum > currentNumInputs
            % make mux larger and have proper number of inports
            set_param(muxSource,'Inputs',num2str(inputNum));
            set_param(muxSource,'Position',[80 38 85 38+20*inputNum]);
            for i = (currentNumInputs+1):inputNum
                % create new number of needed inports
                block_name = [block '/In' num2str(i)];
                add_block('simulink/Sources/In1',block_name);
                % set position of current inport block
                set_param(block_name,'Position',[20 38+(i-1)*30 50 52+(i-1)*30]);
                % connect to proper mux port
                add_line(block,['In' num2str(i) '/1'],['Mux/' num2str(i)])
            end
        elseif inputNum < currentNumInputs
            for i = currentNumInputs:-1:(inputNum+1)
                % delete unnecessary lines
                delete_line(block,['In' num2str(i) '/1'],['Mux/' num2str(i)])
                % delete unnecessary inports
                currentBlock = [block '/In' num2str(i)];
                delete_block(currentBlock);
            end
            % resize mux and change inputs
            set_param(muxSource,'Inputs',num2str(inputNum));
            set_param(muxSource,'Position',[80 38 85 38+20*inputNum]);
        else
            % do nothing. Number of inports equals number of desired
            % inports.
        end
    end
    
    if outputNum ~= -1
        currentNumOutputs = length(portInfo.Outport);
        demuxSource = [block '/Demux'];
        if outputNum > currentNumOutputs
            % make demux larger and have proper number of outports
            set_param(demuxSource,'Outputs',num2str(outputNum));
            set_param(demuxSource,'Position',[280 38 285 38+20*outputNum]);
            for i = (currentNumOutputs+1):outputNum
                % create new number of needed outports
                block_name = [block '/Out' num2str(i)];
                add_block('simulink/Sinks/Out1',block_name);
                % set position of current outport block
                set_param(block_name,'Position',[320 38+(i-1)*30 350 52+(i-1)*30]);
                % connect to proper demux port
                add_line(block,['Demux/' num2str(i)],['Out' num2str(i) '/1'])
            end
        elseif outputNum < currentNumOutputs
            for i = currentNumOutputs:-1:(outputNum+1)
                % delete unnecessary lines
                delete_line(block,['Demux/' num2str(i)],['Out' num2str(i) '/1']);
                % delete unnecessary outports
                currentBlock = [block '/Out' num2str(i)];
                delete_block(currentBlock);
            end
            % resize demux and change outputs
            set_param(demuxSource,'Outputs',num2str(outputNum));
            set_param(demuxSource,'Position',[280 38 285 38+20*outputNum]);
        else
            % do nothing. Number of outports equals number of desired
            % outports.
        end
    end
end