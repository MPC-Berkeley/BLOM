function isquadratic = BLOM_CheckIfQuadratic(P,K)

% Check if the constraints are linear

Pconstr = P(any(K(2:end,:),1),:);
if ~BLOM_CheckIfLinear(Pconstr)
    isquadratic = false;
    return
end

Pcost = P(any(K(1,:),1),:);



%Alowed powers are 0, 1 and 2
illegal_powers =  Pcost(Pcost ~= 0 & Pcost ~= 1 & Pcost ~= 2);
if ~isempty(illegal_powers)
        isquadratic = false;
        return;
end

P_pattern = spones(Pcost);
nnz_P_per_row = sum(P_pattern, 2);
sum_of_powers = sum(abs(Pcost), 2);

if any((nnz_P_per_row > 2) | (sum_of_powers > 2 ) )
    isquadratic = false;
else
    isquadratic = true;
end