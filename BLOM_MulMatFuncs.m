function y = BLOM_MulMatFuncs(f,g,size_f,size_g,same_vars)


m=1;
for k=1:size_g(2) % g columns
    for i=1:size_f(1) % f rows
        if (same_vars)
            Aki = sparse(1,size(f.A,2));
        else
            Aki = sparse(1,size(f.A,2)+size(g.A,2));
        end
        Cki = 0;
        for j=1:size_g(1) % f columns (and g rows)
            [A C] = mul_poly(f.A,g.A,f.C((j-1)*size_f(1)+i,:),g.C((k-1)*size_g(1)+j,:),same_vars);
            [Aki,Cki]= sum_poly(Aki,Cki,A,C);
            
        end
        [Aki,Cki]= BLOM_Filter_poly(Aki,Cki);
        
        y.A{i,k} = Aki;
        y.C{i,k} = Cki;
    end
end

A = [];
for j=1:size(y.A,2)
    for i=1:size(y.A,1)
        A = [A ; y.A{i,j}];
    end
end

k=1;

n=0;
C = [];
for j=1:size(y.A,2)
    for i=1:size(y.A,1)
        C(k,n+[1:size(y.C{i,j},2)]) = y.C{i,j};
        k=k+1;
        n = n + size(y.C{i,j},2);
    end
end

[A,C] = BLOM_Filter_poly(A,C);

y.A =A;
y.C =C;

return 


function  [A , C ] = sum_poly(A1,C1,A2,C2)

    A = [A1 ; A2];
    C = [C1 C2];
%     [A C] = Filter_poly(A,C);
    


function [A C] = mul_poly(A1,A2,C1,C2,same_vars)

c1_idx = find(C1);
c2_idx = find(C2);
if isempty(c1_idx) || isempty(c2_idx)
    if (same_vars)
        A = sparse(1,size(A1,2),0);
    else
        A = sparse(1,size(A1,2)+size(A2,2),0);
    end
    C = 0;
    return;
end

if (same_vars)
    A = sparse(size(A1,1)*size(A2,1),size(A1,2),0);
else
    A = sparse(size(A1,1)*size(A2,1),size(A1,2)+size(A2,2),0);
end
C = zeros(1,size(A1,1)*size(A2,1));
k=1;

for i=c1_idx
    for j=c2_idx
    if (same_vars)
        A(k,:) =   A1(i,:) + A2(j,:);
    else
        A(k,:) =  [ A1(i,:) A2(j,:)];
    end
        C(k)   =  C1(i)*C2(j);
        k=k+1;
    end
end
% [A C] = Filter_poly(A,C);


function [Aunique Cunique] = Filter_poly(A,C)

[Aunique, I, J] = unique(A,'rows');
idxs = 1:size(Aunique,1);
sources = idxs(J);
Cunique = zeros(size(C,1),size(Aunique,1));
for j=1:size(C,1)
    for i=idxs
        Cunique(j,i) = sum(C(j,sources == i));
    end
end