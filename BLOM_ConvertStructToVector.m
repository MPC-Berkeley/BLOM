function vec = BLOM_ConvertStructToVector(all_names,data)


vec = zeros(length(all_names),1);
for i=1:length(all_names)
    
    % for assignment from structure to a vector, any filed out of multiple
    % fields is O.K. So just extract the fist one.
    [name R] = strtok(all_names{i},';');
    
    
    dot_idx = strfind(name,'.') ;
    base_name = name(1:dot_idx(1)-1);
    port_number = str2double(name(dot_idx(1)+4:dot_idx(2)-1));
    time_index = str2double(name(dot_idx(2)+2:end));
    if strcmp(base_name(1:3),'BL_')
        fname = base_name(4:end);
    else
        fname = base_name;
    end
    vec(i) = data.(fname)(time_index ,port_number) ;
    
    
    
end

