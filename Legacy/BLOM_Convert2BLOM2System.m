%> @file Convert2BLOM2system.m
%> @brief converts all BLOM1 blocks within system to BLOM2 blocks.  Does
%>  nothing if no blocks need to be converted.  BLOMLib, Simulink, and
%>  system
%>
%> @param BLOM1system String for name of system.  BLOM1system must be open
%>
%======================================================================

function BLOM_Convert2BLOM2System(BLOM1System)
    load_system('BLOM_Lib') % make sure BLOM2 library is loaded
    subsystems = find_system({BLOM1System});
    subsystems(1) = [];
    for k = 1:length(subsystems)
        BLOM_Convert2BLOM2Block(subsystems{k});
    end
end