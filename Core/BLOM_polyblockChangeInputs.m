%> @file BLOM_polyblockChangeInputs.m
%> @brief callback function for Polyblock that changes the number of inputs
%> to the block based on what the user specified
%>
%> More detailed description of the problem.
%>
%> @param block the current handle of the block
%>
%======================================================================

function BLOM_polyblockChangeInputs(block)
    params = get_param(block,'polyInput');
    paramSplit = regexp(params,',','split');
    for i = 1:length(paramSplit)
        %add inputs for each parameter
    end
    add_block('simulink/Sources/In1',block)
end