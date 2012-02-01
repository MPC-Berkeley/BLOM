function [A C] = BLOM_ConstTrMatConst(const1,X,const2)
% BLOM_ConstMatConst
% Returns matrix of the form Z = c1*X'*c2
% where c1, c2 are constant matrices, X is a variable matrix

f = BLOM_MulConstVar(eye(size(X,1)),X);
[Axt Cxt] = BLOM_TranposeMat(f.A,f.C,size(X));
f.A = Axt;
f.C = Cxt;
c1.A = sparse(1,size(f.A,2),0);
c1.C = const1(:);

c2.A = sparse(1,size(f.A,2),0);
c2.C = const2(:);
y = BLOM_MulMatFuncs(c1,f,[size(const1,1) ,  size(const1,2)],size(X'),true);
z = BLOM_MulMatFuncs(y,c2,[size(const1,1) size(X',2) ],[size(const2,1) ,  size(const2,2)],true);

A = z.A;
C = z.C;

