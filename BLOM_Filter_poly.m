function [Aunique Cunique] = BLOM_Filter_poly(A,C)

[Aunique, I, J] = unique(A,'rows');
idxs = 1:size(Aunique,1);
sources = idxs(J);
Cunique = zeros(size(C,1),size(Aunique,1));
for j=1:size(C,1)
    for i=idxs
        Cunique(j,i) = sum(C(j,sources == i));
    end
end