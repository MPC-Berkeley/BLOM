function str = GetBlockName(name)
if nargin == 0
    name = gcb;
end

idx = strfind(name,'/');
name(idx) = '_';
name = name(idx(1)+1:end);
idx = strfind(name,' ');
name(idx) = '_';

str = name;