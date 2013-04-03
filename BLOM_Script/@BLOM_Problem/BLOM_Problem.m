classdef BLOM_Problem < handle
    properties
        Pt = sparse(0, 0); % store transpose of P internally
        K = sparse(1, 0); % first row of K is cost function
        lb = -inf(0, 1);
        ub = inf(0, 1);
    end
    properties (Dependent = true)
        P
        cost
    end
    
    methods
        function newvar = newVariable(problem, numrows, numcolumns)
            % create a new variable of size numrows * numcolumns
            if nargin < 3
                % assume column vector if numcolumns not given (this is
                % different than behavior of ones(n), zeros(n), sdpvar(n))
                numcolumns = 1;
                if nargin < 2
                    % assume scalar if numrows not given
                    numrows = 1;
                elseif numel(numrows) > 1
                    if numel(numrows) > 2
                        warning('extra entries in numrows ignored')
                    end
                    numcolumns = numrows(2);
                    numrows = numrows(1);
                end
            elseif numel(numrows) > 1 || numel(numcolumns) > 1
                error('numrows and numcolumns should be scalars if both given')
            end
            % add new rows to problem.Pt and save indices in newvar
            idx = size(problem.Pt, 1) + reshape(1:(numrows*numcolumns), ...
                numrows, numcolumns);
            newvar = BLOM_Variable(problem, idx);
            problem.Pt = blkdiag(problem.Pt, sparse(numrows*numcolumns, 0));
            % add entries to lb and ub for new variables
            problem.lb = [problem.lb; -inf(numrows*numcolumns, 1)];
            problem.ub = [problem.ub; inf(numrows*numcolumns, 1)];
        end
        
        function P = get.P(problem)
            P = problem.Pt';
        end
        function set.P(problem, P)
            if size(P, 2) ~= size(problem.Pt, 1)
                error(['P must have the same number of columns as ' ...
                    'the number of rows in problem.Pt'])
            elseif size(P, 1) ~= size(problem.K, 2)
                error(['P must have the same number of rows as ' ...
                    'the number of columns in problem.K'])
            end
            problem.Pt = P';
        end
        
        function cost = get.cost(problem)
            costK = problem.K(1, :);
            costPt = problem.Pt;
            % remove terms that are not used in the cost
            costPt = costPt(:, costK ~= 0);
            costK  = costK(1, costK ~= 0);
            cost = BLOM_Expression(problem, costPt, costK);
        end
        function set.cost(problem, cost)
            if isempty(cost)
                % set cost to zero
                problem.K(1, :) = 0;
                return
            elseif isnumeric(cost)
                % constant cost
                cost = BLOM_Expression(problem, ...
                    sparse(numel(problem.lb), 1), cost(:), false);
            else
                % convert to BLOM_Expression (refreshes size of cost.Pt)
                cost = BLOM_Expression(cost);
            end
            if ~isequal(cost.problem, problem)
                error('cost input does not belong to this BLOM_Problem')
            elseif size(cost.K, 1) > 1
                error('cost must be a scalar')
            else
                problem.Pt = [problem.Pt, cost.Pt, cost.auxPt];
                % clear any existing cost from K
                problem.K(1, :) = 0;
                % add cost.K to end of first row, pad zeros below it
                problem.K = blkdiag([problem.K, blkdiag(cost.K, ...
                    sparse(size(problem.Pt, 1) - 1, 0))], cost.auxK);
            end
        end
        
        function setBounds(problem, expression, lb, ub)
            % modify bounds if expression is a variable,
            % otherwise introduce slack variable and set its bounds
            % (unless lb == 0 and ub == 0, then just modify P and K)
            % maybe output the slack variable?
            if isnumeric(expression)
                error('cannot set bounds for a constant quantity')
            elseif ~isequal(expression.problem, problem)
                error('expression input does not belong to this BLOM_Problem')
            end
            if nargin < 3 || isempty(lb)
                lb = -inf;
            end
            if nargin < 4 || isempty(ub)
                ub = inf;
            end
            if isa(expression, 'BLOM_Variable')
                boundsize = size(expression.idx);
            else
                % convert to BLOM_Expression (refreshes size of expression.Pt)
                expression = BLOM_Expression(expression);
                boundsize = [size(expression.K, 1), 1];
            end
            if numel(lb) == 1
                lb = lb*ones(boundsize);
            end
            if numel(ub) == 1
                ub = ub*ones(boundsize);
            end
            if ~isequal(size(lb), boundsize)
                error('lb must be a scalar or the same size as expression')
            elseif ~isequal(size(ub), boundsize)
                error('ub must be a scalar or the same size as expression')
            end
            if isa(expression, 'BLOM_Variable')
                idx = expression.idx;
                problem.lb(idx) = max(problem.lb(idx), lb);
                problem.ub(idx) = min(problem.ub(idx), ub);
            else
                % add to P & K where lb == ub == 0,
                % otherwise introduce slack vars
                if any(lb ~= 0) || any(ub ~= 0)
                    slack_idx = find((lb ~= 0) | (ub ~= 0));
                    num_slacks = numel(slack_idx);
                    slackvar = problem.newVariable(num_slacks, 1);
                    idx = slackvar.idx;
                    problem.lb(idx) = max(problem.lb(idx), lb);
                    problem.ub(idx) = min(problem.ub(idx), ub);
                    if any(problem.lb(idx) > problem.ub(idx))
                        warning('infeasible bound set')
                    end
                    newPt = blkdiag(expression.Pt, speye(num_slacks));
                    newK = [expression.K, sparse(slack_idx, 1:num_slacks, ...
                        -1, size(expression.K, 1), num_slacks)];
                    expression.auxPt = blkdiag(expression.auxPt, ...
                        sparse(num_slacks, 0));
                else
                    newPt = expression.Pt;
                    newK = expression.K;
                end
                problem.Pt = [problem.Pt, newPt, expression.auxPt];
                problem.K = blkdiag(problem.K, newK, expression.auxK);
            end
        end
        
        function removeUnusedTerms(problem)
            % remove columns of Pt corresponding to empty columns of K
            used_terms = full(any(problem.K ~= 0, 1))';
            problem.Pt = problem.Pt(:, used_terms);
            problem.K  = problem.K(:, used_terms);
        end
        
        function combineDuplicateTerms(problem)
            % find the unique columns of Pt
            problem.removeUnusedTerms;
            
        end
    end
end
