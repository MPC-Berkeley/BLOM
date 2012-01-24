function [y gy ] = GenericCostFun(x,cost,costGrad)

y = cost(x);
gy = costGrad(x);
