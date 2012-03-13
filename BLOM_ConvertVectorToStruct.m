function data = BLOM_ConvertVectorToStruct(all_names,vec)

base_name = cell(1,length(all_names));
port_number = zeros(1,length(all_names));
time_index = zeros(1,length(all_names));

for i=1:length(all_names)
    dot_idx = strfind(all_names{i},'.') ;
    base_name{i} = all_names{i}(1:dot_idx(1)-1);
    port_number(i) = str2double(all_names{i}(dot_idx(1)+4:dot_idx(2)-1));
    time_index(i) = str2double(all_names{i}(dot_idx(2)+2:end));
end

[names, I , J] = unique(base_name);


for i=length(all_names):-1:1 % backward loop prevents reallocating
    if strcmp(names{J(i)}(1:3),'BL_')
        name = names{J(i)}(4:end);
    else
        name = names{J(i)};
    end
    data.(name)(time_index(i) ,port_number(i)) = vec(i);
end