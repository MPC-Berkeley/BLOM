%> @file BLOM_polyblockChangeInputs.m
%> @brief callback function for Polyblock that changes the number of inputs
%> to the block based on user specified P and K matricies
%>
%> More detailed description of the problem.
%>
%> @param block the current handle of the block
%>
%======================================================================

function BLOM_PolyblockChangeInputs(block,P,K,inputs,inputScalar,outputScalar)
    % may want to keep this code to label ports on the block itself
    % paramSplit = regexp(inputs,',','split');
    
    if inputScalar
        if outputScalar
            BLOM_ChangeNumberOfPorts(block,size(P,2),size(K,1));
        else
            BLOM_ChangeNumberOfPorts(block,size(P,2),1);
        end
    else
        if outputScalar
            BLOM_ChangeNumberOfPorts(block,1,size(K,1));
        else
            BLOM_ChangeNumberOfPorts(block,1,1);
        end
    end
end