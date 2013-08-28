%> @file BLOM_SetDataLogging.m
%> @brief sets data logging to on for all variables and labels them
%>
%> @param BLOMSystem String for name of system.  BLOMSystem must be open
%>
%======================================================================

function BLOM_SetDataLogging(BLOMSystem)
    set_param(BLOMSystem, 'StrictBusMsg', 'ErrorLevel1')
	set_param(gcs, 'SignalLoggingSaveFormat', 'ModelDataLogs');

    subsystems = find_system({BLOMSystem});
    systemNameLength = length(BLOMSystem);
    subsystems(1) = [];
    for k = 1:length(subsystems)
        blockType = get_param(subsystems{k}, 'BlockType');
        if ~(strcmp(blockType, 'Mux') || strcmp(blockType, 'Demux') ||  ...
                strcmp(blockType, 'Inport') || strcmp(blockType, 'From') || ...
                (strcmp(blockType, 'SubSystem') && isempty(get_param(subsystems{k}, 'ReferenceBlock'))))

            name = get_param(subsystems{k}, 'Name');
            if ~isvarname(name)
                warning([name ' is not a valid variable name.  Please rename this block and try again.']);
                continue
            end
            
            parent = get_param(subsystems{k}, 'Parent');
            if ~strcmp(parent, BLOMSystem)
                parent(1:systemNameLength+1) = [];
                name = [parent '/' name];
                name(name=='/') = '_';
            end
            portHandles = get_param(subsystems{k}, 'PortHandles');
            outports = portHandles.Outport;

            for j = 1:length(outports)
                set_param(outports(j), 'DataLoggingNameMode', 'custom');
                set_param(outports(j), 'DataLogging', 'on');
                set_param(outports(j), 'DataLoggingName', sprintf('%s_%d', name, j));
            end
        end
    end
end
