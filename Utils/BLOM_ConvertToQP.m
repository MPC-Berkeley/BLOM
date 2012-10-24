function [qp vars_included] = BLOM_ConvertToQP(P, K, lb, ub, options)
% This function converts the following optimization problem:
% min  BLOM_EvalPolyBlock(P, K(1,:), x)
%  x   s.t. lb <= x <= ub,
%      and  BLOM_EvalPolyBlock(P, K(2:end,:), x) = 0 (elementwise)
% into a quadratic program of the form:
% min  0.5*z'*qp.H*z + qp.f'*z
%  z   s.t. qp.lb <= z <= qp.ub,
%      and  qp.Aeq*z = qp.beq,
% where z = x(vars_included). Exits with an error if the original
% optimization problem cannot be converted into a quadratic program.

n = size(P,2); % number of variables in original problem
if length(lb) ~= n
    error('lb is length %d, should be length %d', length(lb), n)
elseif length(ub) ~= n
    error('ub is length %d, should be length %d', length(ub), n)
end

if nargin > 4 && isfield(options, 'vars_not_to_remove')
    vars_not_to_remove = options.vars_not_to_remove;
else
    vars_not_to_remove = false(n,1);
end
%[P K vars_removed] = BLOM_SubstituteIntoCost(P, K, lb, ub, vars_not_to_remove);
vars_removed = false(n,1);

Pvals = nonzeros(P);
if any(Pvals ~= 1 & Pvals ~= 2)
    error(['Problem contains at least one entry in P ' ...
        'not equal to 0, 1, or 2, so is not quadratic.'])
end
term_degrees = sum(P, 2);
if any(term_degrees > 2)
    error(['Problem contains at least one polynomial term (row of P) ' ...
        'with degree greater than 2, so is not quadratic.'])
end
constant_terms   = (term_degrees == 0);
linear_terms     = (term_degrees == 1);
quadratic_terms  = (term_degrees == 2);
constraint_terms = any(K(2:end,:), 1)';
if any(constraint_terms & quadratic_terms)
    error(['At least one constraint is quadratic or bilinear, so cannot ' ...
        'convert to QP. This may be fixable when BLOM_SubstituteIntoCost is done.'])
    %error(['At least one quadratic or bilinear constraint remains after ' ...
    %    'substituting simple constraints into cost, so cannot convert to QP.'])
end
vars_included = ~vars_removed;
n_included = sum(vars_included);
Ptrans = P';
Ptrans_linear = Ptrans(:, linear_terms);
Ptrans_quadratic = Ptrans(:, quadratic_terms);
qp.lb = lb(vars_included);
qp.ub = ub(vars_included);
qp.Aineq = sparse(0, n_included);
qp.bineq = zeros(0, 1);
qp.Aeq = K(2:end, linear_terms) * Ptrans_linear';
qp.beq = -sum(K(2:end, constant_terms), 2);
qp.f = Ptrans_linear * K(1, linear_terms)';

quadratic_cost_coeffs = K(1, quadratic_terms)';
[var term pow] = find(Ptrans_quadratic);
var_bilinear = var(pow == 1);
var_quadratic = var(pow == 2);
term_bilinear = term(pow == 1);
term_quadratic = term(pow == 2);
if ~isequal(term_bilinear(1:2:end), term_bilinear(2:2:end))
    error('Unexpected ordering in bilinear terms in P')
end
var1 = [var_quadratic; var_bilinear(1:2:end)];
var2 = [var_quadratic; var_bilinear(2:2:end)];
vals = quadratic_cost_coeffs([term_quadratic; term_bilinear(1:2:end)]);
qp.H = sparse(var1, var2, vals, n_included, n_included);
qp.H = qp.H + qp.H'; % add symmetric triangular part (and double diagonals)

if nargin < 5 || ~isfield(options, 'check_convexity') || options.check_convexity
    [L D perm] = ldl(qp.H, 'vector');
    mineig = min(eig(D));
    if mineig < 0
        warning(['Cost matrix is indefinite, min(eig(D)) from LDL ' ...
            'factorization is %g'], mineig)
    end
end
