function [A, C] = BLOM_ConstMatConst(const1,X,const2)
% BLOM_ConstMatConst
% Returns matrix of the form Z = c1*X*c2
% where c1, c2 are constant matrices, X is a variable matrix

f = BLOM_MulConstVar(const1,X);
p.A = sparse(1,size(f.A,2),0);
p.C = const2(:);
y = BLOM_MulMatFuncs(f,p,[size(const1,1) ,  size(X,2)],[size(const2,1) size(const2,2)],true);

A = y.A;
C = y.C;

