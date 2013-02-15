function [H, LambdaMap] = BLOM_EvalHessian(P, K, x, Lambda)
% This function evaluates the Hessian of the PolyBlock function represented
% by the exponent matrix P and coefficient matrix K, at the column vector x
% with Lagrange multipliers Lambda. The Hessian is a square symmetric
% matrix with dimension equal to the number of columns of P, and only the
% lower triangular part is returned by this function. The first row of K
% represents the cost function. Lambda should be a column vector of length
% equal to the number of rows in K, or the number of rows in K minus 1 (in
% which case the multiplier on the cost is taken equal to 1).
% If Lambda is not given, then a vector of all 1's is assumed. If x is not
% given then this function returns the sparsity pattern of the Hessian.
% The optional second output argument, LambdaMap, indicates which
% constraints contribute to which nonzeros in the Hessian.
% LambdaMap is of dimension m by nnz, where m is the number of rows of
% K and nnz is the number of nonzero elements in the Hessian.

P = P(any(K,1),:); % remove unused terms
K = K(:,any(K,1));

n = size(P,2); % number of variables
r = size(P,1); % number of (used) polynomial terms
m = size(K,1); % number of constraints plus one

if nargin < 4
    Lambda = ones(m, 1); % default value
elseif size(Lambda, 2) > 1
    Lambda = Lambda'; % column vector by convention
end
if length(Lambda) < m-1 || length(Lambda) > m
    error('Lambda is length %d, should be length %d', length(Lambda), m-1)
elseif length(Lambda) == m-1
    Lambda = [1; Lambda]; % prepend a 1 for the cost row
end
if nargin > 2 && length(x) ~= n
    error('x is length %d, should be length %d', length(x), n)
end

[Pcols Prows Pvals] = find(P'); % transpose so nonzeros are in row-major order
P_pattern = spones(P);
nnz_P_per_row = sum(P_pattern, 2);
nnz_P_prev_rows = [0; cumsum(full(nnz_P_per_row))];
nonlinear_terms = find(any(P ~= P_pattern, 2) | (nnz_P_per_row > 1))';

if nargin > 2
    % evaluate polyblock function at given point x
    expbool  = (Pvals == BLOM_FunctionCode('exp'));
    logbool  = (Pvals == BLOM_FunctionCode('log'));
    sinbool  = (Pvals == BLOM_FunctionCode('sin'));
    cosbool  = (Pvals == BLOM_FunctionCode('cos'));
    tanhbool = (Pvals == BLOM_FunctionCode('tanh'));
    atanbool = (Pvals == BLOM_FunctionCode('atan'));
    erfbool  = (Pvals == BLOM_FunctionCode('erf'));
    
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
    
    vx2deriv = (Pvals - 1).*Pvals.*(x(Pcols).^(Pvals - 2)); % 2nd derivatives of powers
    vx2deriv(Pvals == 1) = 0; % fix NaNs from linear powers at x = 0
    vx2deriv(expbool)  = vx(expbool); % 2nd derivatives of exponentials
    vx2deriv(logbool)  = -x(Pcols(logbool)).^(-2); % 2nd derivatives of logarithms
    vx2deriv(sinbool)  = -vx(sinbool); % 2nd derivatives of sines
    vx2deriv(cosbool)  = -vx(cosbool); % 2nd derivatives of cosines
    vx2deriv(tanhbool) = -2*vx(tanhbool).*vxderiv(tanhbool); % 2nd derivatives of hyperbolic tangents
    vx2deriv(atanbool) = -2*x(Pcols(atanbool)).*vxderiv(atanbool).^2; % 2nd derivatives of arctangents
    vx2deriv(erfbool)  = -2*x(Pcols(erfbool)).*vxderiv(erfbool); % 2nd derivatives of error functions
end

% The Hessian pattern of each term is the outer product of the vector of
% nonzeros in that term with itself, except that the diagonal is only
% nonzero if the exponent for that variable is not equal to 1
H_pattern = spones(tril(P_pattern'*P_pattern)); % just lower triangular part
H_diag = any(P ~= P_pattern, 1)';
H_pattern = spdiags(H_diag, 0, H_pattern);
if nargin == 2 && nargout == 1
    % just return Hessian sparsity pattern
    H = H_pattern;
    return
end

[Hrows Hcols] = find(H_pattern);
nnz_H_per_col = sum(H_pattern, 1)';
nnz_H_prev_cols = [0; cumsum(full(nnz_H_per_col))];
nnz_H_allterms = sum(nnz_P_per_row.*(nnz_P_per_row+1)/2) - nnz(P == 1);

PtransPattern = spones(P)'; % save transpose because getting a
% column of a sparse matrix is much faster than getting a row

% Hmap contains the nonzero values of each term's Hessian located according to
% the overall Hessian sparsity pattern, and is of dimension nnz(H_pattern) by r
Hmaprows = zeros(nnz_H_allterms, 1);
Hmapcols = zeros(nnz_H_allterms, 1);
Hmapvals = zeros(nnz_H_allterms, 1);
Hmapend = 0;
% put the elements of Hmap in the correct places now that H_pattern is known
if nargin > 2
    for i=nonlinear_terms
        term_nzinds = (nnz_P_prev_rows(i)+1 : nnz_P_prev_rows(i+1));
        term_vars = Pcols(term_nzinds);
        vx_term = vx(term_nzinds);
        vxderiv_term = vxderiv(term_nzinds);
        vx2deriv_term = vx2deriv(term_nzinds);
        
        termHcolumn = PtransPattern(:, i); % preallocate for this term
        for j=1:length(term_vars)
            % columns of this term's Hessian
            colj = term_vars(j);
            factors = vx_term;
            factors(j) = vxderiv_term(j); % derivative of this column
            nextderiv = vxderiv_term;
            nextderiv(j) = vx2deriv_term(j); % 2nd derivative on diagonal
            for k=j:length(term_vars)
                % rows of this term's Hessian
                if nextderiv(k) == 0
                    termHcolumn(term_vars(k)) = 0;
                else
                    elem_factors = factors;
                    elem_factors(k) = nextderiv(k); % mixed 2nd derivative
                    termHcolumn(term_vars(k)) = prod(elem_factors);
                end
            end
            termHcolreduced = termHcolumn(H_pattern(:, colj));
            if any(termHcolreduced)
                inds = find(termHcolreduced);
                Hmaprange = Hmapend + 1 : Hmapend + length(inds);
                Hmaprows(Hmaprange) = nnz_H_prev_cols(colj) + inds;
                Hmapcols(Hmaprange) = i;
                Hmapvals(Hmaprange) = termHcolreduced(inds);
                Hmapend = Hmapend + length(inds);
            end
        end
    end
    Hmap = sparse(Hmaprows(1:Hmapend), Hmapcols(1:Hmapend), ...
        Hmapvals(1:Hmapend), nnz(H_pattern), r);
    Hvals = Hmap * (K' * Lambda);
    H = sparse(Hrows, Hcols, Hvals, n, n);
    if nargout > 1
        LambdaMap = K * Hmap';
    end
else
    % no x to evaluate at, just return sparsity patterns of H and LambdaMap
    H = H_pattern;
    termHdiagonals = (P' ~= PtransPattern);
    for i=nonlinear_terms
        term_nzinds = (nnz_P_prev_rows(i)+1 : nnz_P_prev_rows(i+1));
        term_vars = Pcols(term_nzinds);
        
        termHcolumn = PtransPattern(:, i);
        for j=1:length(term_vars)
            % columns of this term's Hessian
            colj = term_vars(j);
            if H_diag(colj)
                prev = termHcolumn(colj);
                termHcolumn(colj) = termHdiagonals(colj, i);
                termHcolreduced = termHcolumn(H_pattern(:, colj));
                termHcolumn(colj) = prev; % restore entry for next column
            else
                termHcolreduced = termHcolumn(H_pattern(:, colj));
            end            
            if any(termHcolreduced)
                inds = find(termHcolreduced);
                Hmaprange = Hmapend + 1 : Hmapend + length(inds);
                Hmaprows(Hmaprange) = nnz_H_prev_cols(colj) + inds;
                Hmapcols(Hmaprange) = i;
                Hmapend = Hmapend + length(inds);
            end
        end
    end
    Hmap = sparse(Hmaprows, Hmapcols, 1, nnz(H_pattern), r);
    LambdaMap = spones(spones(K) * Hmap');
end
