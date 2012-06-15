function [A, C] = BLOM_MatMat(A_X,C_X,size_x,A_Y,C_Y,size_y,same_var)
% BLOM_MatMat
% Returns matrix of the form Z = X*Y
% where  X, Y are given as a matrix functions
% A_X,C_X - A and C matrices of X function
% A_Y,C_Y - A and C matrices of Y function
% size_x,size_y - sizes of X and Y matrices respectively
% same_var - if true, X and Y are functions of the same variables, in this
%            case, the result A matrix has the same number of columns as X
%            and Y.
%            if false, X and Y are functions of different variables. In
%            this case the result A matrix has number of columns as a sum
%            of X and Y.


x.A = A_X;
x.C = C_X;
y.A = A_Y;
y.C = C_Y;

y = BLOM_MulMatFuncs(x,y,size_x,size_y,same_var);

A = y.A;
C = y.C;

