function y = BLOM_EvalPolyBlock(P, K, x)
% This function evaluates the PolyBlock function represented by the
% exponent matrix P and coefficient matrix K, at the column vector x

P = P(any(K,1),:); % remove unused terms
K = K(:,any(K,1));

if 1 % sparse version
    [Prows Pcols Pvals] = find(P);
    vx = x(Pcols).^Pvals; % powers of input variables
    expbool = (Pvals == BLOM_FunctionCode('exp'));
    vx(expbool) = exp(x(Pcols(expbool))); % exponentials
    logbool = (Pvals == BLOM_FunctionCode('log'));
    vx(logbool) = log(x(Pcols(logbool))); % logarithms
    prods = ones(size(P,1),1); % preallocate products vector
    for v = 1:length(Pvals)
        prods(Prows(v)) = prods(Prows(v)) * vx(v); % compute products for each term
    end
else % dense version
    xrep = ones(size(P,1),1)*x'; % transpose x and repeat for each term
    vx = xrep.^P; % powers of input variables
    expbool = (P == BLOM_FunctionCode('exp'));
    vx(expbool) = exp(xrep(expbool)); % exponentials
    logbool = (P == BLOM_FunctionCode('log'));
    vx(logbool) = log(xrep(logbool)); % logarithms
    prods = prod(vx, 2); % compute products for each term
end
y = K*prods; % output vector
