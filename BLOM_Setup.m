function BLOM_Setup


BLOM_dir = which('BLOM_setup');
BLOM_dir = BLOM_dir(1:strfind(BLOM_dir,'BLOM_Setup.m')-1);


switch (computer)
    case {'PCWIN','PCWIN64'}
        build_ipopt = false;        
    case  {'GLNX86','GLNXA64'}
        build_ipopt = true; 
        makefile = 'Makefile';
    case {'MACI','MACI64'}
        build_ipopt = true;
        makefile = 'Makefile.mac';
end



if (build_ipopt)
    ipopt_dir = uigetdir('','Pick an IPOPT folder that holds Lib directory, press cancel for no IPOPT');
    if (ipopt_dir ~= 0)
        cur_dir = pwd;
        cd(BLOM_dir);
        cd('BLOM_Ipopt');
        eval(['! make -f ' makefile ' all IPOPTPATH=' ipopt_dir]);
        
        if (~exist('BLOM_NLP','file'))
            warning('Compilation of BLOM_NLP failed. Check the screen for errors');
        else
            disp('---------------------------------');
            disp('BLOM_NLP was succesfully compiled');
            disp('---------------------------------');
        end
        cd(cur_dir)
        
        
    end
end