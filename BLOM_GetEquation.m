function [equation] = BLOM_GetEquation(ModelSpec,type, num)


all_names = cell(length(ModelSpec.all_names));
for i=1:length(ModelSpec.all_names)
    name = strtok(ModelSpec.all_names{i},';');
    all_names{i} = name(4:end);
end

for k=1:length(num)

switch(type)
    case 'eq'
        count = 0;
        for i=1:length(ModelSpec.AAs)
            
            count = count  + size(ModelSpec.Cs{i},1);
            if count >= num(k)
                break;
            end
            
        end        
        equation{k} = CreatePolyFromMatrix(ModelSpec.AAs{i},ModelSpec.Cs{i}(end-(count-num(k)),:),all_names,'low');
    case 'ineq'
        count = 0;
        for i=1:length(ModelSpec.ineq.AAs)
            
            count = count  + size(ModelSpec.ineq.Cs,1);
            if count >= num(k)
                break;
            end
            
        end        
        equation{k} = CreatePolyFromMatrix(ModelSpec.ineq.AAs{i},ModelSpec.ineq.Cs{i}(end-(count-num(k)),:),all_names,'low'); 
    case 'cost'
        equation{k} = CreatePolyFromMatrix(ModelSpec.cost.A,ModelSpec.cost.C,all_names,'low'); 
        
end
        

end