%> @file Convert2BLOM2Block.m
%> @brief converts a single BLOM1 block into BLOM2 block.  Does nothing if
%>  block does not need to be converted.  BLOMLib, Simulink, and file 
%>  containing block must all be open 
%>
%> @param BLOM1Block block. can be form of name or handle
%>
%======================================================================
function BLOM_Convert2BLOM2Block(BLOM1Block)

    refBlock = get_param(BLOM1Block, 'ReferenceBlock');
    position = get_param(BLOM1Block, 'Position');
    parent = get_param(BLOM1Block, 'Parent');
    name = get_param(BLOM1Block, 'Name');
    orientation = get_param(BLOM1Block, 'Orientation');
    block = [parent '/' name];

    switch refBlock
        case 'MPCMdlLib/UserFriendly PolyBlock'
            scalarInput = strcmp(get_param(BLOM1Block, 'input_type'), 'Scalar');
            scalarOutput = strcmp(get_param(BLOM1Block, 'output_type'), 'Scalar');
            P_f = get_param(BLOM1Block, 'A');
            K_f = get_param(BLOM1Block, 'C');
            P_g = get_param(BLOM1Block, 'gA');
            K_g = get_param(BLOM1Block, 'gC');
            inportNames = get_param(BLOM1Block, 'port_names');

            delete_block(BLOM1Block);
            add_block('BLOM_Lib/GeneralPolyblock', block);            
                     
            set_param(block, 'P_f', P_f);
            set_param(block, 'K_f', K_f);
            set_param(block, 'P_g', P_g);
            set_param(block, 'K_g', K_g);
            set_param(block, 'inputs', inportNames);
            if scalarInput
                set_param(block, 'inputScalar', 'on');
            else
                set_param(block, 'inputScalar', 'off');
            end
            if scalarOutput
                set_param(block, 'outputScalar', 'on');
            else 
                set_param(block, 'outputScalar', 'off');
            end
            
            set_param(block, 'Orientation', orientation);
            set_param(block, 'Position', position);

            
        case 'MPCMdlLib/Constraint'
            positive = strcmp(get_param(BLOM1Block, 'CompareSign'), '> 0');
            
            
            add_block('BLOM_Lib/Bound', block);

            set_param(block, 'initial_step', 'on');
            set_param(block, 'final_step', 'on');
            set_param(block, 'intermediate_step', 'on');
            if positive
                set_param(block, 'ub', 'inf');
                set_param(block, 'lb', '0');
            else
                set_param(block, 'ub', '0');
                set_param(block, 'lb', '-inf');
            end
            
            delete_block(BLOM1Block);
            set_param(block, 'Orientation', orientation);
            set_param(block, 'Position', position);
            
        case 'MPCMdlLib/state'
            initialCondition = get_param(BLOM1Block, 'init_val');
            
            delete_block(BLOM1Block);
            add_block('simulink/Discrete/Unit Delay', block);
            set_param(block, 'Position', position);
            
            set_param(block, 'X0', initialCondition);
        
        case 'MPCMdlLib/Continuous state'
            initialCondition = get_param(BLOM1Block, 'init_val');
            
            delete_block(BLOM1Block);
            add_block('simulink/Continuous/Integrator', block);
            set_param(block, 'Orientation', orientation);
            set_param(block, 'Position', position);
            
            set_param(block, 'InitialCondition', initialCondition);
            
        case 'MPCMdlLib/Cost functon'
            delete_block(BLOM1Block);
            add_block('BLOM_Lib/DiscreteCost', block);
            set_param(block, 'Orientation', orientation);
            set_param(block, 'Position', position);

        case 'MPCMdlLib/Input var'
            delete_block(BLOM1Block);
            add_block('BLOM_Lib/InputFromSimulink', block);
            set_param(block, 'Position', position);
            
        case 'MPCMdlLib/External Var'
            delete_block(BLOM1Block);
            add_block('BLOM_Lib/ExternalFromSimulink', block);
            set_param(block, 'Orientation', orientation);
            set_param(block, 'Position', position);
            
        case 'MPCMdlLib/PolyBlock'
            P = eval(get_param(BLOM1Block, 'A'));
            K = eval(get_param(BLOM1Block, 'C'));
            numOutputs = size(K,1);
            
            inputCols = P(:,1:end-numOutputs);
            outputCols = P(:,end-numOutputs+1:end);
            outputIdx = any(outputCols == 1, 2);
            if any(any(outputCols ~= 1 & outputCols ~= 0))
                error('Invalid values in output columns of A.  Must be only 1 or 0')
            end
            if any(any(inputCols ~= 0,2) & any(outputCols ~= 0,2))
                error('Output variables should only be identity.  Should have no relation to inputs');
            end
            if any(sum(outputCols,2) > 1)
                error('Invalid outputs.  Rows in A should not relate multiple outputs.')
            end
            if any(sum(outputCols,1) > 1)
                error('Multiple rows of A for same output variable is redundant. Please remove such that there is only one instance.')
            end
            if any(any(K(:,outputIdx) ~= -1 & K(:,outputIdx) ~= 0))
                error('All entries in C corresponding to output variables in A should be -1');
            end
            P(outputIdx,:) = [];
            P(:,end-numOutputs+1:end) = [];
            K(:,outputIdx) = [];

            scalarInput = strcmp(get_param(BLOM1Block, 'input_type'), 'Scalar');
            scalarOutput = strcmp(get_param(BLOM1Block, 'output_type'), 'Scalar');
            
            delete_block(BLOM1Block);
            add_block('BLOM_Lib/Polyblock', block);            
            set_param(block, 'P', P);
            set_param(block, 'K', K);
            if scalarInput
                set_param(block, 'inputScalar', 'on');
            else
                set_param(block, 'inputScalar', 'off');
            end
            if scalarOutput
                set_param(block, 'outputScalar', 'on');
            else 
                set_param(block, 'outputScalar', 'off');
            end
            set_param(block, 'Orientation', orientation);
            set_param(block, 'Position', position);

    end

end