function [c ceq gradc gradceq ] = GenericConstrFun(x,neq,eq,gneq,geq)

c = neq(x);
ceq = eq(x);
gradc = gneq(x);
gradceq = geq(x);
