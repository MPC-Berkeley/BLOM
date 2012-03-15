function str = PrintProblem(all_names, AAs ,  Cs , ineq,cost)

str = {};

s = 'J = ';
s = [s  CreatePolyFromMatrix(cost.A,cost.C,all_names)];
str{end+1} = s ;



for i=1:length(AAs)
    for j=1:size(Cs{i},1)
    str{end+1} = [ CreatePolyFromMatrix(AAs{i},Cs{i}(j,:),all_names) ' = 0'] ;
    end
end

for i = 1:length(ineq.AAs)
    for j=1:size(ineq.Cs{i},1)
        str{end+1} = [ CreatePolyFromMatrix(ineq.AAs{i},ineq.Cs{i}(j,:),all_names) ' < 0'] ;
    end
end