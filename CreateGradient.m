function grad = CreateGradient(A,C)
grad.AAs = {};
grad.Cs = {};

for i = 1:size(A,2) % for all variables
    grad.AAs{i} = sparse(size(A,1),size(A,2),0);
    grad.Cs{i} = sparse(size(C),1,0);
    for j=1:size(A,1) % for all terms
        if A(j,i) == 0
            continue;
        end
        
        if ~isinf(A(j,i))  
            grad.AAs{i}(j,:) = A(j,:);
            grad.AAs{i}(j,i) = grad.AAs{i}(j,i)-1;
            grad.Cs{i}(j) = C(j)*A(j,i);
        else
            if (A(j,i) > 0 ) % exp
                grad.AAs{i}(j,:) = A(j,:);
                grad.Cs{i}(j) = C(j);
            else  % log
                grad.AAs{i}(j,:) = A(j,:);
                grad.AAs{i}(j,i) = -1;
                grad.Cs{i}(j) = C(j);
            end

        end
    end
end