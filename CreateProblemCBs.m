function [cost_name costGrad_name eqconstr_name eqconstrGrad_name neqconstr_name neqconstrGrad_name ...
          Aeq Beq A B ]...
    =CreateProblemCBs(name,all_names, AAs ,  Cs , ineq,cost,separate_linear)



for i = 1:length(all_names)
    idx_names{i} = sprintf('x(%d)',i);
end


str =  CreatePolyFromMatrix(cost.A,cost.C,idx_names);



cost_name = [name 'Cost'];
PrintFunction(cost_name,str);

str = {};

grad = CreateGradient(cost.A,cost.C);
str = {};
for i = 1:length(grad.AAs)
        str{end+1} =  CreatePolyFromMatrix(grad.AAs{i},grad.Cs{i},idx_names);
end


costGrad_name = [name 'CostGrad'];
PrintFunction(costGrad_name,str);


if (separate_linear)
 [ Aeq , Beq ,   AAs , Cs ] = ExtractLinearPart(AAs,Cs);
else
    Aeq = [];
    Beq = [];
end

str = {};
for i=1:length(AAs)
    for j=1:size(Cs{i},1)
        str{end+1} = CreatePolyFromMatrix(AAs{i},Cs{i}(j,:),idx_names)  ;
    end
end



eqconstr_name = [name 'Eq'];
PrintFunction(eqconstr_name,str);

clear grad;
for i=1:length(AAs)
    grad{i} = CreateGradient(AAs{i},Cs{i});
end

for i=1:length(AAs)
    for j=1:length(idx_names)
        str{i,j} = CreatePolyFromMatrix(grad{i}.AAs{j},grad{i}.Cs{j},idx_names);
    end
end

eqconstrGrad_name = [name 'EqGrad'];
PrintFunction(eqconstrGrad_name,str);

if (separate_linear)
 [ A , B ,   ineq.AAs , ineq.Cs ] = ExtractLinearPart(ineq.AAs,ineq.Cs);
else
    A = [];
    B = [];
end

str = {};
for i = 1:length(ineq.AAs)
        for j=1:size(ineq.Cs{i},1)
            str{end+1} =  CreatePolyFromMatrix(ineq.AAs{i},ineq.Cs{i}(j,:),idx_names)  ;
        end
end

neqconstr_name = [name 'NEq'];
PrintFunction(neqconstr_name,str)


clear grad;
for i=1:length(ineq.AAs)
    grad{i} = CreateGradient(ineq.AAs{i},ineq.Cs{i});
end

for i=1:length(ineq.AAs)
    for j=1:length(idx_names)
        str{i,j} = CreatePolyFromMatrix(grad{i}.AAs{j},grad{i}.Cs{j},idx_names);
    end
end

neqconstrGrad_name = [name 'NEqGrad'];
PrintFunction(neqconstrGrad_name,str)


function PrintFunction(name,str)

f = fopen([name '.m'],'wt');

fprintf(f,'%%#eml\n');
fprintf(f,'function y = %s(x)\n\n',name);

if ~iscell(str)
    str = {str};
end

MinimumMatrixSizeForSpare = 40;

if size(str,1)*size(str,2) == 1  
    % do not initialize
elseif size(str,1)*size(str,2) > MinimumMatrixSizeForSpare 
    % sparse
    fprintf(f,'y = sparse(%d,%d,0);\n\n',size(str,2),size(str,1));
else
    fprintf(f,'y = zeros(%d,%d);\n\n',size(str,2),size(str,1));
end

for i=1:size(str,1)
    for j=1:size(str,2)
        if ~strcmp(str{i,j},'0')
            fprintf(f,'y(%d,%d) = %s ; \n',j,i,str{i,j});
        end
    end
end

fclose(f);


function [ Matrix , B,   nlAAs , nlCs ] = ExtractLinearPart(AAs,Cs)

l=1;
nl=1;
Matrix = [];
B = [];
nlAAs = {};
nlCs  = {};
for i = 1:length(AAs)
    if IsLinear(AAs{i})
        for j=1:size(Cs{i},1)
            Matrix(l,:) = Cs{i}(j,:)*AAs{i};
            idx = find(sum(AAs{i}')==0);% find constant terms
            if isempty(idx)
                B(l) = 0;
            else
                B(l) = -sum(Cs{i}(j,idx));
            end

            l=l+1;
        end
       else
           nlAAs{nl} = AAs{i};
           nlCs{nl} = Cs{i};
           nl = nl+1;
       end
           
    
end

Matrix = sparse(Matrix);

function is_linear = IsLinear(A)
if (max(sum(A'))<=1)&& (min(sum(A')) > -1 )
    is_linear = true;
else
    is_linear = false;
end
    
 