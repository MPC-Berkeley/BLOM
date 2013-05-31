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
            
            set_param(block, 'Position', position);
            
        case 'MPCMdlLib/Constraint'
            positive = strcmp(get_param(BLOM1Block, 'CompareSign'), '> 0');
            
            delete_block(BLOM1Block);
            add_block('BLOM_Lib/Bound', block);
            set_param(block, 'Position', position);

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
            
        case 'MPCMdlLib/state'
            initialCondition = get_param(BLOM1Block, 'init_val');
            
            delete_block(BLOM1Block);
            add_block('simulink/Discrete/Unit Delay', block);
            set_param(block, 'Position', position);
            
            set_param(block, 'X0', initialCondition);
            
        case 'MPCMdlLib/Cost functon'
            delete_block(BLOM1Block);
            add_block('BLOM_Lib/DiscreteCost', block);
            set_param(block, 'Position', position);
            

    end


end