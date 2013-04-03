classdef (InferiorClasses = {?BLOM_Variable}) BLOM_Expression
    properties
        problem % handle of parent problem
        Pt % store transpose of P internally - make hidden eventually?
        K
        specialFunction
        auxPt % P transpose for any auxiliary equalities introduced by expression
        auxK % K for any auxiliary equalities introduced by expression
    end
    properties (Dependent = true)
        P
    end
    
    methods
        function expr = BLOM_Expression(problem, Pt, K, specialFunction)
            if isa(problem, 'BLOM_Problem')
                expr.problem = problem;
                if nargin < 2
                    expr.Pt = sparse(size(problem.Pt, 1), 0);
                elseif size(Pt, 1) ~= size(problem.Pt, 1)
                    error('Pt must have the same number of rows as problem.Pt')
                else
                    expr.Pt = Pt;
                end
                if nargin < 3
                    expr.K = sparse(0, size(expr.Pt, 2));
                elseif size(K, 2) ~= size(expr.Pt, 2)
                    error('Pt and K must have the same number of columns')
                else
                    expr.K = K;
                end
                if nargin < 4
                    expr.specialFunction = any(ismember(nonzeros(expr.Pt), ...
                        BLOM_FunctionCode('all_codes')));
                else
                    expr.specialFunction = specialFunction;
                end
                expr.auxPt = sparse(size(problem.Pt, 1), 0);
                expr.auxK = sparse(0, 0);
            elseif isa(problem, 'BLOM_Expression')
                expr = problem;
                % update the number of rows in expr.Pt if is not equal
                % to size(expr.problem.Pt, 1), meaning variables have
                % been added to the parent BLOM_Problem object since
                % the expression was first created
                numvars1 = size(expr.Pt, 1);
                numvars2 = size(expr.problem.Pt, 1);
                if numvars1 > numvars2
                    error('expression somehow has more rows in Pt than its parent problem')
                elseif numvars1 < numvars2
                    expr.Pt = blkdiag(expr.Pt, sparse(numvars2 - numvars1, 0));
                end
                numvars3 = size(expr.auxPt, 1);
                if numvars3 > numvars2
                    error('expression somehow has more rows in auxPt than its parent problem')
                elseif numvars3 < numvars2
                    expr.auxPt = blkdiag(expr.auxPt, sparse(numvars2 - numvars3, 0));
                end
            elseif isa(problem, 'BLOM_Variable')
                var = problem;
                idx = var.idx(:);
                expr = BLOM_Expression(var.problem, sparse(idx, ...
                    1:numel(idx), 1, numel(var.problem.lb), ...
                    numel(idx)), speye(numel(idx)), false);
            elseif isa(problem, 'BLOM_Constraint')
                error('invalid operation on BLOM_Constraint object')
            else
                error('first input is not a BLOM object')
            end
        end
        
        function P = get.P(expr)
            P = expr.Pt';
        end
        function expr = set.P(expr, P)
            if size(P, 2) ~= size(expr.problem.Pt, 1)
                error(['P must have the same number of columns as ' ...
                    'the number of rows in expr.problem.Pt'])
            elseif size(P, 1) ~= size(expr.K, 2)
                error(['P must have the same number of rows as ' ...
                    'the number of columns in expr.K'])
            end
            expr.Pt = P';
            expr.specialFunction = any(ismember(nonzeros(expr.Pt), ...
                BLOM_FunctionCode('all_codes')));
        end
        
        function expr = removeUnusedTerms(expr)
            % remove columns of Pt corresponding to empty columns of K
            used_terms = full(any(expr.K ~= 0, 1))';
            if ~any(used_terms) && ~isempty(expr.K) && ~isempty(expr.Pt)
                % no terms used, output one column of zero terms
                expr.K  = 0*expr.K(:, 1);
                expr.Pt = 0*expr.Pt(:, 1);
            else
                expr.Pt = expr.Pt(:, used_terms);
                expr.K  = expr.K(:, used_terms);
            end
            expr.specialFunction = any(ismember(nonzeros(expr.Pt), ...
                BLOM_FunctionCode('all_codes')));
        end
        
        function expr = combineDuplicateTerms(expr)
            % find the unique columns of Pt
            expr = expr.removeUnusedTerms;
            
        end
        
        % plus, multiplication, division, and powers are in separate files
        
        function out = uplus(in1) % +x
            out = in1;
        end
        
        function out = minus(in1, in2)
            out = plus(in1, uminus(in2));
        end
        
        function out = uminus(in1) % -x
            out = in1;
            out.K = -out.K;
        end
        
        function out = mrdivide(in1, in2)
            if ~isnumeric(in1)
                in1 = BLOM_Expression(in1);
            end
            if isnumeric(in2)
                size2 = size(in2);
            else
                in2 = BLOM_Expression(in2);
                size2 = [size(in2.K, 1), 1];
            end
            if max(size2) > 1
                error(['matrix division not supported, use ./ for ' ...
                    'elementwise division'])
            end
            out = in1*(in2^(-1));
        end
        
        function out = rdivide(in1, in2)
            if ~isnumeric(in1)
                in1 = BLOM_Expression(in1);
            end
            if ~isnumeric(in2)
                in2 = BLOM_Expression(in2);
            end
            out = in1.*(in2.^(-1));
        end
        
        function out = ctranspose(in1)
            warning(['all BLOM_Expression objects are considered vectors, ' ...
                'ctranspose does nothing'])
            out = in1;
        end
        
        function out = transpose(in1)
            warning(['all BLOM_Expression objects are considered vectors, ' ...
                'transpose does nothing'])
            out = in1;
        end
        
        function out = sqrt(in1)
            out = power(in1, 0.5);
        end
        
        function out = exp(in1)
            out = power(in1, BLOM_FunctionCode('exp'));
        end
        
        function out = log(in1)
            out = power(in1, BLOM_FunctionCode('log'));
        end
        
        function out = sin(in1)
            out = power(in1, BLOM_FunctionCode('sin'));
        end
        
        function out = cos(in1)
            out = power(in1, BLOM_FunctionCode('cos'));
        end
        
        function out = tanh(in1)
            out = power(in1, BLOM_FunctionCode('tanh'));
        end
        
        function out = atan(in1)
            out = power(in1, BLOM_FunctionCode('atan'));
        end
        
        function out = erf(in1)
            out = power(in1, BLOM_FunctionCode('erf'));
        end
        
        % dot
        % prod
        % sum
        % horzcat
        % vertcat
        % repmat
        % reshape?
        % subsref
        % subsasgn
        % le (output BLOM_Constraint object)
        % ge (output BLOM_Constraint object)
        % eq (output BLOM_Constraint object)
        
        % other trig, etc functions either derived or implicitly?
    end
end
