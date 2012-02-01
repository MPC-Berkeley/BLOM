function [A, C] = A_BLOM_MatMul(const1,X,const2,Y,const3,same_vars)
% A_BLOM_MatMul_c1Xc2Yc3
% Returns matrix of the form Z = c1*X*c2*Y*c3
% where c1, c2, c3 are constant matrices, X and Y are variable matrices

if nargin<6
    same_vars = false;
end

f = MulConstVar(const1,X);
g = MulConstVar(const2,Y);
w = BLOM_MulMatFuncs(f,g,[size(const1,1) ,  size(X,2)],[size(const2,1) size(Y,2)],same_vars);
p.A = sparse(1,size(w.A,2),0);
p.C = const3(:);
y = BLOM_MulMatFuncs(w,p,[size(const1,1) ,  size(Y,2)],[size(const3,1) size(const3,2)],true);

A = y.A;
C = y.C;



function f = MulConstVar(const,X)

f.A = speye(size(X,1)*size(X,2));
f.C = sparse(size(X,2)*size(const,1),size(X,1)*size(X,2));

for j=1:size(X,2)
    for i=1:size(const,1)
        f.C(i+(j-1)*size(const,1),(1:size(const,2))+(j-1)*size(const,2))= const(i,:);
    end
end
