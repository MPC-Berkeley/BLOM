function str = PrintProblem(all_names, AAs ,  Cs , ineq,cost)

str = {};

s = 'J = ';
s = [s  CreatePolyFromMatrix(cost.A,cost.C,all_names)];
str{end+1} = s ;



for i=1:length(AAs)
    str{end+1} = [ CreatePolyFromMatrix(AAs{i},Cs{i},all_names) ' = 0'] ;
end

for i = 1:length(ineq.AAs)
        str{end+1} = [ CreatePolyFromMatrix(ineq.AAs{i},ineq.Cs{i},all_names) ' > 0'] ;
end