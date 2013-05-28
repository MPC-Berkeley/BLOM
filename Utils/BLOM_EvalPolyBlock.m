function y = BLOM_EvalPolyBlock(P, K, x)
% This function evaluates the PolyBlock function represented by the
% exponent matrix P and coefficient matrix K, at the column vector x

persistent BLOM_FunctionCodes
if isempty(BLOM_FunctionCodes)
    BLOM_FunctionCodes = BLOM_FunctionCode('codes_struct');
end

P = P(any(K,1),:); % remove unused terms
K = K(:,any(K,1));

if 1 % sparse version
    [Prows Pcols Pvals] = find(P);
    expbool  = (Pvals == BLOM_FunctionCodes.exp);
    logbool  = (Pvals == BLOM_FunctionCodes.log);
    sinbool  = (Pvals == BLOM_FunctionCodes.sin);
    cosbool  = (Pvals == BLOM_FunctionCodes.cos);
    tanhbool = (Pvals == BLOM_FunctionCodes.tanh);
    atanbool = (Pvals == BLOM_FunctionCodes.atan);
    erfbool  = (Pvals == BLOM_FunctionCodes.erf);
    
    vx = x(Pcols).^Pvals; % powers of input variables
    vx(expbool)  = exp(x(Pcols(expbool))); % exponentials
    vx(logbool)  = log(x(Pcols(logbool))); % logarithms
    vx(sinbool)  = sin(x(Pcols(sinbool))); % sines
    vx(cosbool)  = cos(x(Pcols(cosbool))); % cosines
    vx(tanhbool) = tanh(x(Pcols(tanhbool))); % hyperbolic tangents
    vx(atanbool) = atan(x(Pcols(atanbool))); % arctangents
    vx(erfbool)  = erf(x(Pcols(erfbool))); % error functions
    
    prods = ones(size(P,1),1); % preallocate products vector
    for v = 1:length(Pvals)
        prods(Prows(v)) = prods(Prows(v)) * vx(v); % compute products for each term
    end
else % dense version
    xrep = ones(size(P,1),1)*x'; % transpose x and repeat for each term
    expbool  = (P == BLOM_FunctionCodes.exp);
    logbool  = (P == BLOM_FunctionCodes.log);
    sinbool  = (P == BLOM_FunctionCodes.sin);
    cosbool  = (P == BLOM_FunctionCodes.cos);
    tanhbool = (P == BLOM_FunctionCodes.tanh);
    atanbool = (P == BLOM_FunctionCodes.atan);
    erfbool  = (P == BLOM_FunctionCodes.erf);
    
    vx = xrep.^P; % powers of input variables
    vx(expbool)  = exp(xrep(expbool)); % exponentials
    vx(logbool)  = log(xrep(logbool)); % logarithms
    vx(sinbool)  = sin(xrep(sinbool)); % sines
    vx(cosbool)  = cos(xrep(cosbool)); % cosines
    vx(tanhbool) = tanh(xrep(tanhbool)); % hyperbolic tangents
    vx(atanbool) = atan(xrep(atanbool)); % arctangents
    vx(erfbool)  = erf(xrep(erfbool)); % error functions
    
    prods = prod(vx, 2); % compute products for each term
end
y = K*prods; % output vector
