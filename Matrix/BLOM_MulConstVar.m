function f = BLOM_MulConstVar(const,X)

f.A = speye(size(X,1)*size(X,2));
f.C = sparse(size(X,2)*size(const,1),size(X,1)*size(X,2));

for j=1:size(X,2)
    for i=1:size(const,1)
        f.C(i+(j-1)*size(const,1),(1:size(const,2))+(j-1)*size(const,2))= const(i,:);
    end
end
