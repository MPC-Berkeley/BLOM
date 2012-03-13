function vec = BLOM_ConvertStructToVector(all_names,data)

% base_name = cell(1,length(all_names));
% port_number = zeros(1,length(all_names));
% time_index = zeros(1,length(all_names));
vec = zeros(length(all_names),1);
for i=1:length(all_names)
    dot_idx = strfind(all_names{i},'.') ;
    base_name = all_names{i}(1:dot_idx(1)-1);
    port_number = str2double(all_names{i}(dot_idx(1)+4:dot_idx(2)-1));
    time_index = str2double(all_names{i}(dot_idx(2)+2));
    if strcmp(base_name(1:3),'BL_')
        name = base_name(4:end);
    else
        name = base_name;
    end
    vec(i) = data.(name)(time_index ,port_number) ;
end

