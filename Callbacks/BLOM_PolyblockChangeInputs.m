%> @file BLOM_polyblockChangeInputs.m
%> @brief callback function for Polyblock that changes the number of inputs
%> to the block based on user specified P and K matricies
%>
%> More detailed description of the problem.
%>
%> @param block the current handle of the block
%>
%======================================================================

function BLOM_PolyblockChangeInputs(block)
    % may want to keep this code to label ports on the block itself
    params = get_param(block,'inputs');
    paramSplit = regexp(params,',','split');
    
    P = eval(get_param(block,'P'));
    K = eval(get_param(block,'K'));
    input_as_scalar = get_param(block,'inputScalar');
    output_as_scalar = get_param(block,'outputScalar');
    if input_as_scalar
        if output_as_scalar
            BLOM_ChangeNumberOfPorts(block,size(P,2)-size(K,1),size(K,1));
        else
            BLOM_ChangeNumberOfPorts(block,size(P,2)-size(K,1),1);
        end
    else
        if output_as_scalar
            BLOM_ChangeNumberOfPorts(block,1,size(K,1));
        else
            BLOM_ChangeNumberOfPorts(block,1,1);
        end
    end
    
    %% following code just added for testing purposes. to be removed later.
    constantName = [block '/Constant'];
    K = eval(get_param(block,'K'));
    valueString = ['ones(' num2str(length(K(:,1))) ',1)'];
    set_param(constantName,'Value',valueString);
end