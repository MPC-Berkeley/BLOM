function [P K vars_removed] = ...
    BLOM_SubstituteIntoCost(P, K, lb, ub, vars_not_to_remove)
% This function finds and eliminates unbounded variables that appear in
% the cost function and only one simple constraint, except for variables
% specified in the optional last input. The last output vars_removed is
% a logical array with length equal to the original number of variables.

P = P(any(K,1),:); % remove unused terms
K = K(:,any(K,1));
% combine duplicate rows of P
[P idx1 idx2] = unique(P, 'rows');
perm2 = sparse(1:length(idx2), idx2, 1, length(idx2), size(P,1));
K = K*perm2;

n = size(P,2); % number of variables in original problem
if length(lb) ~= n
    error('lb is length %d, should be length %d', length(lb), n)
elseif length(ub) ~= n
    error('ub is length %d, should be length %d', length(ub), n)
end
if nargin < 5
    vars_not_to_remove = false(n,1);
elseif length(vars_not_to_remove) ~= n || max(vars_not_to_remove) > 1
    % assume this was given as a list of indices
    vars_not_to_remove = sparse(vars_not_to_remove, 1, true, n, 1);
end
vars_removed = ~any(P,1)';
if any(vars_removed)
    msg = sprintf('Original problem has %d unused variables.', sum(vars_removed));
    if nargout < 3
        warning(msg)
    else
        warning([msg ' vars_removed output includes these.'])
    end
    if any(vars_removed & vars_not_to_remove)
        warning('At least one of the indicated vars_not_to_remove is unused!')
    end
end
% don't remove bounded variables
vars_not_to_remove = vars_not_to_remove | (lb ~= -inf) | (ub ~= inf);

vars_tried = vars_removed; % initialize
P_pattern = spones(P);
J_count = (spones(K) * P_pattern)';
vars_in_cost = J_count(:,1);
J_count = J_count(:,2:end)';
constraints_per_var = sum(spones(J_count), 1)';
while 1
    var_to_try = find(~vars_tried & ~vars_not_to_remove & ...
        vars_in_cost & (constraints_per_var == 1), 1);
    if isempty(var_to_try)
        break
    end
    vars_tried(var_to_try) = true;
    [constraint_to_try one count_terms] = find(J_count(:, var_to_try));
    if count_terms > 1
        % more than one term in this constraint depends on this variable
        continue
    end
    % more work to do here...
end
P(:,vars_removed) = [];

