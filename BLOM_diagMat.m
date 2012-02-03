function [A C] = BLOM_diagMat(Ain,Cin,size_x)

A = Ain;

idx = 1:size_x(1)*size_x(2);

mat_idx = reshape(idx,size_x(1),size_x(2));
mat_idx = diag(mat_idx);
C = Cin(mat_idx(:),:);
