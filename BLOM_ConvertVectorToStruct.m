function data = BLOM_ConvertVectorToStruct(all_names,vec)

base_name = cell(1,length(all_names));
port_number = zeros(1,length(all_names));
time_index = zeros(1,length(all_names));

k=1;
for i=1:length(all_names)
    [name R] = strtok(all_names{i},';');
    while ~isempty(name)
        dot_idx = strfind(name,'.') ;
        base_name{k} = name(1:dot_idx(1)-1);
        port_number(k) = str2double(name(dot_idx(1)+4:dot_idx(2)-1));
        time_index(k) = str2double(name(dot_idx(2)+2:end));
        vec_idx(k) = i;
        k =k+1;
        [name R] = strtok(R,';');
    end
end

[names, I , J] = unique(base_name);


for i=length(base_name):-1:1 % backward loop prevents reallocating
    if strcmp(names{J(i)}(1:3),'BL_')
        name = names{J(i)}(4:end);
    else
        name = names{J(i)};
    end
    data.(name)(time_index(i) ,port_number(i)) = vec(vec_idx(i));
end