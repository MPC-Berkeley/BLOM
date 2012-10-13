function y = BLOM_EvalPolyBlock(P, K, x)
% This function evaluates the PolyBlock function represented by the
% exponent matrix P and coefficient matrix K, at the column vector x

P = P(any(K,1),:); % remove unused terms
K = K(:,any(K,1));

if 1 % sparse version
    [Prows Pcols Pvals] = find(P);
    expbool  = (Pvals == BLOM_FunctionCode('exp'));
    logbool  = (Pvals == BLOM_FunctionCode('log'));
    sinbool  = (Pvals == BLOM_FunctionCode('sin'));
    cosbool  = (Pvals == BLOM_FunctionCode('cos'));
    tanhbool = (Pvals == BLOM_FunctionCode('tanh'));
    
    vx = x(Pcols).^Pvals; % powers of input variables
    vx(expbool)  = exp(x(Pcols(expbool))); % exponentials
    vx(logbool)  = log(x(Pcols(logbool))); % logarithms
    vx(sinbool)  = sin(x(Pcols(sinbool))); % sines
    vx(cosbool)  = cos(x(Pcols(cosbool))); % cosines
    vx(tanhbool) = tanh(x(Pcols(tanhbool))); % hyperbolic tangents
    
    prods = ones(size(P,1),1); % preallocate products vector
    for v = 1:length(Pvals)
        prods(Prows(v)) = prods(Prows(v)) * vx(v); % compute products for each term
    end
else % dense version
    xrep = ones(size(P,1),1)*x'; % transpose x and repeat for each term
    expbool  = (P == BLOM_FunctionCode('exp'));
    logbool  = (P == BLOM_FunctionCode('log'));
    sinbool  = (P == BLOM_FunctionCode('sin'));
    cosbool  = (P == BLOM_FunctionCode('cos'));
    tanhbool = (P == BLOM_FunctionCode('tanh'));
    
    vx = xrep.^P; % powers of input variables
    vx(expbool)  = exp(xrep(expbool)); % exponentials
    vx(logbool)  = log(xrep(logbool)); % logarithms
    vx(sinbool)  = sin(xrep(sinbool)); % sines
    vx(cosbool)  = cos(xrep(cosbool)); % cosines
    vx(tanhbool) = tanh(xrep(tanhbool)); % hyperbolic tangents
    
    prods = prod(vx, 2); % compute products for each term
end
y = K*prods; % output vector
