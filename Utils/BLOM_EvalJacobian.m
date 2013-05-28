function J = BLOM_EvalJacobian(P, K, x)
% This function evaluates the Jacobian of the PolyBlock function represented
% by the exponent matrix P and coefficient matrix K, at the column vector x.
% The Jacobian has as many rows as K, and as many columns as P. If the
% first row of K represents the cost function, then the first row of the
% output Jacobian represents the gradient of the cost function.

persistent BLOM_FunctionCodes
if isempty(BLOM_FunctionCodes)
    BLOM_FunctionCodes = BLOM_FunctionCode('codes_struct');
end

P = P(any(K,1),:); % remove unused terms
K = K(:,any(K,1));

n = size(P,2); % number of variables
r = size(P,1); % number of (used) polynomial terms

[Pcols Prows Pvals] = find(P'); % transpose so nonzeros are in row-major order
P_pattern = spones(P);
nnz_P_per_row = sum(P_pattern, 2);
nnz_P_prev_rows = [0; cumsum(full(nnz_P_per_row))];
nonlinear_terms = find(any(P ~= P_pattern, 2) | (nnz_P_per_row > 1))';

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

vxderiv = Pvals.*(x(Pcols).^(Pvals - 1)); % derivatives of powers
vxderiv(expbool)  = vx(expbool); % derivatives of exponentials
vxderiv(logbool)  = 1./x(Pcols(logbool)); % derivatives of logarithms
vxderiv(sinbool)  = cos(x(Pcols(sinbool))); % derivatives of sines
vxderiv(cosbool)  = -sin(x(Pcols(cosbool))); % derivatives of cosines
vxderiv(tanhbool) = sech(x(Pcols(tanhbool))).^2; % derivatives of hyperbolic tangents
vxderiv(atanbool) = 1./(x(Pcols(atanbool)).^2 + 1); % derivatives of arctangents
vxderiv(erfbool)  = 2./(sqrt(pi)*exp(x(Pcols(erfbool)).^2)); % derivatives of error functions

% construct Jacobian of vx product vector (same sparsity pattern as P)
prodJacvals = ones(size(Pvals)); % for linear terms, this is 1
for i=nonlinear_terms
    j_range = (nnz_P_prev_rows(i) + 1 : nnz_P_prev_rows(i + 1));
    vxrow = vx(j_range);
    for j=j_range
        factors = vxrow;
        % replace value for variable j with its derivative
        factors(j - nnz_P_prev_rows(i)) = vxderiv(j);
        prodJacvals(j) = prod(factors);
    end
end

J = K * sparse(Prows, Pcols, prodJacvals, r, n);

%{
% verification code using naive method
prodJacvals1 = prodJacvals;
for i=1:length(prodJacvals)
    row = Prows(i);
    row_start = nnz_P_prev_rows(row) + 1;
    row_end = nnz_P_prev_rows(row + 1);
    prodJacvals(i) = prod(vx(row_start:i-1)) * vxderiv(i) * prod(vx(i+1:row_end));
end
if ~isequal(prodJacvals, prodJacvals1)
    error('mismatch')
end
%}
