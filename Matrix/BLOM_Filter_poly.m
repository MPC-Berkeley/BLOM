function [Aunique Cunique] = BLOM_Filter_poly(A,C)

[Aunique, I, J] = unique(A,'rows');
perm = sparse(1:length(J), J, 1, length(J), size(Aunique,1));
Cunique = C*perm;
%{
idxs = 1:size(Aunique,1);
sources = idxs(J);
Cunique_old = zeros(size(C,1),size(Aunique,1));
for j=1:size(C,1)
    for i=idxs
        Cunique_old(j,i) = sum(C(j,sources == i));
    end
end
if ~isequal(Cunique,Cunique_old)
    error('mismatch')
end
%}
