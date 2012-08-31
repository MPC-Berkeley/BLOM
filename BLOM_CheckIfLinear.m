function islinear = BLOM_CheckIfLinear(P)

P_pattern = spones(P);
nnz_P_per_row = sum(P_pattern, 2);
if any(any(P ~= P_pattern, 2) | (nnz_P_per_row > 1))
    islinear = false;
else
    islinear = true;
end