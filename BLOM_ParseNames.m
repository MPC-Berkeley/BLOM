function [port_names_str] = BLOM_ParseNames(port_names,n)

k=0;
fields = textscan(port_names,'%s','Delimiter',',');
for i=1:length(fields{1})
    ranges = textscan(fields{1}{i},'%s%d%d','Delimiter','[]-,');
    if ~isempty(ranges{2})

        for j=ranges{2}:ranges{3}
           str = [ranges{1}{1} num2str(j)];
           k = k+1;
           port_names_str{k} = str;        
           if k>=n 
               break
           end
        end
    else
           k = k+1;
           port_names_str{k} = fields{1}{i};        
    end
    if k>=n
        break
    end
end

if (length(port_names_str) < n)
    for i= length(port_names_str)+1: n
        port_names_str{i} = ['x' sprintf('%d',i)];
    end
end

