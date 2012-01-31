function [A, C] = A_BLOM_MatMul(const1,X,const2,Y,const3)
% A_BLOM_MatMul_c1Xc2Yc3
% Returns matrix of the form Z = c1*X*c2*Y*c3
% where c1, c2, c3 are constant matrices, X and Y are variable matrices

f = MulConstVar(const1,X);
g = MulConstVar(const2,Y);
w = Mul2Funcs(f,g,[size(const1,1) ,  size(X,2)],[size(const2,1) size(Y,2)]);


A = w.A;
C = w.C;



function  w = Mul2Funcs(f,g,size_f,size_g)

m=1;
for k=1:size_g(2) % g columns
    for i=1:size_f(1) % f rows
        Aki = sparse(1,size(f.A,2)+size(g.A,2));
        Cki = 0;
        for j=1:size_g(1) % f columns (and g rows)
            [A C] = mul_poly(f.A,g.A,f.C((j-1)*size_f(1)+i,:),g.C((k-1)+j,:));
            [Aki,Cki]= sum_poly(Aki,Cki,A,C);
            
        end
        [Aki,Cki]= Filter_poly(Aki,Cki);
        
        w.A{i,k} = Aki;
        w.C{i,k} = Cki;
    end
end

A = [];
for i=1:size(w.A,1)
    for j=1:size(w.A,2)
        A = [A ; w.A{i,j}];
    end
end

k=1;

n=0;
C = [];
for i=1:size(w.A,1)
    for j=1:size(w.A,2)
        C(k,n+[1:size(w.C{i,j},2)]) = w.C{i,j};
        k=k+1;
        n = n + size(w.C{i,j},2);
    end
end

[A,C] = Filter_poly(A,C);

w.A =A;
w.C =C;

return 


function  [A , C ] = sum_poly(A1,C1,A2,C2)

    A = [A1 ; A2];
    C = [C1 C2];
%     [A C] = Filter_poly(A,C);
    


function [A C] = mul_poly(A1,A2,C1,C2)

c1_idx = find(C1);
c2_idx = find(C2);
if isempty(c1_idx) || isempty(c2_idx)
    A = sparse(1,size(A1,2)+size(A2,2),0);
    C = 0;
    return;
end
A = sparse(size(A1,1)*size(A2,1),size(A1,2)+size(A2,2),0);
C = zeros(1,size(A1,1)*size(A2,1));
k=1;

for i=c1_idx
    for j=c2_idx
        A(k,:) =  [ A1(i,:) A2(j,:)];
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

function f = MulConstVar(const,X)

f.A = speye(size(X,1)*size(X,2));
f.C = sparse(size(X,2)*size(const,1),size(X,1)*size(X,2));

for j=1:size(X,2)
    for i=1:size(const,1)
        f.C(i+(j-1)*size(const,1),(1:size(const,2))+(j-1)*size(const,2))= const(i,:);
    end
end
