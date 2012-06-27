function vec = BLOM_ConvertStructToVector(all_names,data)

if isstruct(all_names)
    % can also input all_names_struct with precomputed vectorization info
    all_fields = all_names.all_fields;
    vec_idx = all_names.vec_idx;
    vec = nan(length(all_names.terms_so_far)-1,1);
else
    % number of ';' is number of multiple names
    num_terms = cellfun(@length, strfind(all_names,';')) + 1;
    terms_so_far = [0, cumsum(num_terms)];
    all_fields = textscan([all_names{:}],'BL_%sOut%dt%d','Delimiter','.;');
    vec_idx = zeros(terms_so_far(end),1); % preallocate vec_idx
    vec_idx(terms_so_far(1:end-1)+1) = 1:length(all_names); % first of each
    twoterms = find(num_terms == 2);
    vec_idx(terms_so_far(twoterms)+2) = twoterms; % 2nd of each
    multiterms = find(num_terms > 2); % multiples, should be fewer of these
    for i = 1:length(multiterms)
        vec_idx(terms_so_far(multiterms(i))+2 : ...
            terms_so_far(multiterms(i)+1)) = multiterms(i);
    end
    vec = nan(length(all_names),1);
end
fnames = fieldnames(data);
for i=1:length(fnames)
    % find matching signals for this field name, but only for the time and
    % port numbers that are saved in the data struct
    matches = find(strcmp(fnames{i}, all_fields{1}) & ...
        all_fields{2} <= size(data.(fnames{i}),2) & ...
        all_fields{3} <= size(data.(fnames{i}),1));
    if isempty(matches)
        error(['Could not find signal named ' fnames{i} ' in all_names'])
    end
    % only use first signal for each match
    matches = matches([true; diff(vec_idx(matches)) ~= 0]);
    % for each field name, convert times and port #'s into 1d indices
    inds_i = sub2ind(size(data.(fnames{i})), all_fields{3}(matches), ...
        all_fields{2}(matches));
    vec(vec_idx(matches)) = data.(fnames{i})(inds_i);
end

%{
vec1 = nan(length(all_names),1);
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
        vec1(i) = data.(fname)(time_index ,port_number) ;
    end
%     end    
end
if ~isequalwithequalnans(vec, vec1)
    disp('mismatch')
end
%}
