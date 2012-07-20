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

function [P,K] = BLOM_convert2PolyBlock(blockHandle)
    % FIX: some way to figure out what type of block it is
    blockType = get_param(blockHandle,'BlockType');
    
    %% switch for different cases. hopefully able to compare numbers
    %% instead of strings
    switch blockType
        
        %% add/subtract 
        case {'Sum','Add'}
            % for n scalar inputs, P = eye(n), K = ones(1,n). for each ith
            % input that is subtracted, that index in K should equal -1
            % instead.
            P = 1;
            K = 1;
            
            
        %% absolute value (only for costs)    
        case 'Abs'
            
            
        %% multiply/divide     
        case {'Product','Divide'}
            % for n scalar inputs, P = ones(1,n), K = [1]. for each ith
            % input to be subtracted, that index in P should be equal -1
            
        %% constant    
        case 'Constant'
            
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