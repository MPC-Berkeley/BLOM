function [HessianPattern, LambdaPattern] = BLOM_HessianPattern(P, K)
% This function calculates and returns the sparsity pattern of the Hessian
% of the PolyBlock function represented by the exponent matrix P and
% coefficient matrix K. The optional second output argument, LambdaPattern,
% indicates which constraints contribute to which nonzeros in the Hessian.
% LambdaPattern is of dimension m by nnz, where m is the number of rows of
% K and nnz is the number of nonzero elements in HessianPattern.

P = P(any(K,1),:); % remove unused terms
K = K(:,any(K,1));

n = size(P,2); % number of variables
r = size(P,1); % number of (used) polynomial terms

% The Hessian pattern of each term is the outer product of the vector of
% nonzeros in that term with itself, except that the diagonal is only
% nonzero if the exponent for that variable is not equal to 1
P_pattern = spones(P);
HessianPattern = spones(tril(P_pattern'*P_pattern)); % just lower triangular part
HessianPattern = spdiags(any(P ~= P_pattern, 1)', 0, HessianPattern);

if nargout > 1
    nnz_Hessian_percolumn = sum(HessianPattern, 1)';
    nnz_Hessian_prevcolumns = [0; cumsum(full(nnz_Hessian_percolumn))];
    
    vars_in_term = sum(P_pattern, 2);
    nnz_Hessian_allterms = sum(vars_in_term.*(vars_in_term+1)/2) - nnz(P == 1);
    
    nonlinear_terms = find(any(P ~= P_pattern, 2) | (vars_in_term > 1));
    nonlinearPtrans = P(nonlinear_terms, :)'; % transpose because getting
    % a column of a sparse matrix is much faster than getting a row
    nonlinearPtransPattern = spones(nonlinearPtrans);
    termHessianDiagonals = nonlinearPtrans ~= nonlinearPtransPattern;
    
    % vxPattern indicates which terms contribute to which nonzeros in the
    % Hessian, and is of dimension nnz by r
    vxPattern = spalloc(nnz(HessianPattern), r, nnz_Hessian_allterms);
    % put the elements of vxPattern in the correct places now that HessianPattern is known
    for i=1:length(nonlinear_terms)
        termNonzeroLocations = nonlinearPtransPattern(:, i);
        termHessianDiagonal = termHessianDiagonals(:, i);
        for j=find(nonlinearPtrans(:, i))'
            termHessianColumn = termNonzeroLocations;
            termHessianColumn(j) = termHessianDiagonal(j);
            vxPattern(nnz_Hessian_prevcolumns(j)+1 : nnz_Hessian_prevcolumns(j+1), ...
                nonlinear_terms(i)) = termHessianColumn(HessianPattern(:, j));
        end
    end
    LambdaPattern = spones(spones(K) * vxPattern');
end
