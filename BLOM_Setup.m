function BLOM_Setup(ipopt_dir)
% Function to easily compile BLOM_NLP for the first time. If called with no
% input arguments, use uigetdir (this will fail in a -nodisplay instance of
% Matlab). Otherwise take first input argument as path to Ipopt.

% mfilename('fullpath') returns the entire path to this script
% the first output of fileparts(ans) gives the path string, same as dirname in unix
BLOM_dir = fileparts(mfilename('fullpath'));

if nargin==0
    disp('Pick an IPOPT folder that holds lib directory, press cancel for no IPOPT')
    ipopt_dir = uigetdir('','Pick an IPOPT folder that holds lib directory, press cancel for no IPOPT');
end
if ~isequal(ipopt_dir, 0) && ~isempty(ipopt_dir)
    if ispc
        % Replace backslashes in BLOM and Ipopt paths with forward slashes
        BLOM_dir = strrep(BLOM_dir, '\', '/');
        ipopt_dir = strrep(ipopt_dir, '\', '/');
        which_shell = questdlg('MinGW or Cygwin?','','MinGW','Cygwin','MinGW');
        if strcmp(which_shell, 'MinGW')
            shell_dir = 'MinGW\msys\1.0';
        elseif strcmp(which_shell, 'Cygwin')
            shell_dir = 'cygwin';
            % Have to translate ipopt_dir into Cygwin /cygdrive/... format
            ipopt_dir = ['`cygpath "' ipopt_dir '"`'];
        else
            error('Invalid selection of shell')
        end
        disp(['Identify ' shell_dir '\bin folder that contains sh.exe'])
        shell_dir = uigetdir(['C:\' shell_dir '\bin'], ...
            ['Identify ' shell_dir '\bin folder that contains sh.exe']);
        if ~exist([shell_dir '\sh.exe'], 'file')
            warning(['Invalid selection of shell_dir=' shell_dir])
        else
            system([shell_dir '\sh --login -c "cd ' BLOM_dir '/BLOM_Ipopt; ' ...
                ipopt_dir '/Ipopt/config.status --file=''Makefile Sparse++/makefile.def''; ' ...
                'make clean; make all"']);
            
            if (~exist('BLOM_NLP.exe','file'))
                warning('Compilation of BLOM_NLP.exe failed. Check the screen for errors');
            else
                disp('-------------------------------------');
                disp('BLOM_NLP.exe was succesfully compiled');
                disp('-------------------------------------');
            end
        end
    else
        cur_dir = pwd;
        cd([BLOM_dir '/BLOM_Ipopt']);
        system([ipopt_dir '/Ipopt/config.status --file="Makefile Sparse++/makefile.def"'])
        system('make clean; make all');
        
        if (~exist('BLOM_NLP','file'))
            warning('Compilation of BLOM_NLP failed. Check the screen for errors');
        else
            disp('---------------------------------');
            disp('BLOM_NLP was succesfully compiled');
            disp('---------------------------------');
        end
        %copyfile([matlabroot '/bin/' computer('arch') '/libmwma57.*'],'.')
        %movefile(ls('libmwma57.*'),strrep(ls('libmwma57.*'),'libmwma57','libhsl'))
        cd(cur_dir)
    end
else
    if ~ischar(ipopt_dir)
        ipopt_dir = num2str(ipopt_dir);
    end
    warning(['Invalid selection of ipopt_dir=' ipopt_dir])
end
