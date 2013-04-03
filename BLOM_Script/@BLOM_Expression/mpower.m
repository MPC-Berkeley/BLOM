function out = mpower(in1, in2)

if ~isnumeric(in2)
    error('variables as exponents is not supported')
elseif numel(in2) > 1
    error('non-scalar powers should use .^')
elseif size(in1.K, 1) > 1
    error('powers of vector variables should use .^')
end
specialFunction = any(ismember(in2, BLOM_FunctionCode('all_codes')));
if nnz(in1.K) <= 1 && ~specialFunction && ~in1.specialFunction
    % single-term expression and no special functions, so can just
    % multiply all exponents in Pt and take powers of K
    out = in1.removeUnusedTerms;
    out.Pt = in2*out.Pt;
    out.K = (out.K)^in2;
else
    % have to introduce new variable to take power of a multi-term
    % expression or special function, or special function of any expression
    newvar = in1.problem.newVariable;
    out = mpower(newvar, in2);
    % set aux Pt and K in output expression for newvar == in1 constraint
    aux = in1 - newvar;
    out.auxPt = [aux.Pt, aux.auxPt];
    out.auxK = blkdiag(aux.K, aux.auxK);
end
