function vec = BLOM_ConvertStructToVector(all_names,data)


vec = nan*ones(length(all_names),1);
for i=1:length(all_names)
    
    % for assignment from structure to a vector, any field out of multiple
    % fields is O.K. So just extract the fist one.
%     [name R] = strtok(all_names{i},';');
%     
%     
%     dot_idx = strfind(name,'.') ;
%     base_name = name(1:dot_idx(1)-1);
%     port_number = str2double(name(dot_idx(1)+4:dot_idx(2)-1));
%     time_index = str2double(name(dot_idx(2)+2:end));
    
    fields = textscan(all_names{i},'BL_%sOut%dt%d','Delimiter','.;');
    
    idx = find(isfield(data,fields{1}),1); % look for the first existing field
    if (isempty(idx))
        continue;
    end
    base_name = fields{1}{idx};
    port_number = fields{2}(idx);
    time_index = fields{3}(idx);
%     if strcmp(base_name(1:3),'BL_')
%         fname = base_name(4:end);
%     else
    fname = base_name;
%     end
%     if (isfield(data,fname))
    if size(data.(fname),1) >= time_index  ...
            && size(data.(fname),2) >= port_number
        vec(i) = data.(fname)(time_index ,port_number) ;
    end
%     end
    
    
end

