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
    params = get_param(block,'polyInput');
    paramSplit = regexp(params,',','split');
    
    P = eval(get_param(block,'polyP'));
    K = eval(get_param(block,'polyK'));
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
end