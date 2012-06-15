function [A C] = BLOM_TraceMat(Ain,Cin,size_x)

A = Ain;

idx = 1:size_x(1)*size_x(2);

mat_idx = reshape(idx,size_x(1),size_x(2));
tr = diag(mat_idx);

C = zeros(1,size(A,1));
for i=1:length(tr)
    C = C + Cin(tr(i),:);
end

[A,C] = BLOM_Filter_poly(A,C);

