function [equation] = BLOM_WhereUsed(ModelSpec,var)


all_names = cell(length(ModelSpec.all_names));
for i=1:length(ModelSpec.all_names)
    name = strtok(ModelSpec.all_names{i},';');
    all_names{i} = name(4:end);
end


C_rows = find(any(ModelSpec.C(:,ModelSpec.A(:,var) ~= 0)'));
equation = {};
for k=C_rows

            
     eq = CreatePolyFromMatrix(ModelSpec.A,ModelSpec.C(k,:),all_names,'low');
     if k == 1
         equation{end+1}   = ['J= ' eq ];
     elseif k >= ModelSpec.ineq_start_C && k <= ModelSpec.ineq_end_C
        equation{end+1}   = [eq '<= 0'];
     elseif k >= ModelSpec.eq_start_C 
        equation{end+1}   = [eq '== 0'];
     end
end
        

