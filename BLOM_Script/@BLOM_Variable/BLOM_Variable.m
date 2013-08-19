% The BLOM_Variable class represents optimization variables according to
% their parent BLOM_Problem object and the index within that problem.
% It should NOT be used directly to create new variables! 
% Use the newVariable method of the BLOM_Problem class to do that.

classdef BLOM_Variable < handle
    properties
        problem % handle of parent problem
        idx % index within parent problem
    end
    properties (Dependent = true)
        value % numerical value (set to nan if no values known)
    end
    
    methods
        function var = BLOM_Variable(problem, idx)
            var.problem = problem;
            var.idx = idx;
        end
        
        function value = get.value(var)
            value = var.problem.x(var.idx);
        end
        function set.value(var, value)
            var.problem.x(var.idx) = value;
        end
        
        function out = subsref(var, sub)
            %{
            if strcmp(sub.type, '.')
                out = var.(sub.subs);
            else
                out = BLOM_Variable(var.problem, subsref(var.idx, sub));
            end
            %}
            out = var;
            for i=1:length(sub)
                if strcmp(sub(i).type, '.')
                    out = out.(sub(i).subs);
                elseif isa(out, 'BLOM_Variable')
                    out = BLOM_Variable(out.problem, subsref(out.idx, sub(i)));
                else
                    out = subsref(out, sub(i));
                end
            end
        end
        
        function out = uplus(in1) % +x
            out = in1;
        end
        
        % most other overloaded operators are defined by converting
        % to a BLOM_Expression object, but have to check for
        % operator(variable, numeric) or operator(numeric, variable)
        function out = plus(in1, in2)
            if isnumeric(in1)
                in1 = BLOM_Expression(in2.problem, ...
                    sparse(numel(in2.problem.x), 1), in1(:), false);
            elseif isnumeric(in2)
                in2 = BLOM_Expression(in1.problem, ...
                    sparse(numel(in1.problem.x), 1), in2(:), false);
            end
            out = plus(BLOM_Expression(in1), BLOM_Expression(in2));
        end
        
        function out = minus(in1, in2)
            if isnumeric(in1)
                in1 = BLOM_Expression(in2.problem, ...
                    sparse(numel(in2.problem.x), 1), in1(:), false);
            elseif isnumeric(in2)
                in2 = BLOM_Expression(in1.problem, ...
                    sparse(numel(in1.problem.x), 1), in2(:), false);
            end
            out = minus(BLOM_Expression(in1), BLOM_Expression(in2));
        end
        
        function out = uminus(in1) % -x
            out = uminus(BLOM_Expression(in1));
        end
        
        function out = mtimes(in1, in2)
            if ~isnumeric(in1)
                in1 = BLOM_Expression(in1);
            end
            if ~isnumeric(in2)
                in2 = BLOM_Expression(in2);
            end
            out = mtimes(in1, in2);
        end
        
        function out = times(in1, in2)
            if ~isnumeric(in1)
                in1 = BLOM_Expression(in1);
            end
            if ~isnumeric(in2)
                in2 = BLOM_Expression(in2);
            end
            out = times(in1, in2);
        end
        
        function out = mrdivide(in1, in2)
            if ~isnumeric(in1)
                in1 = BLOM_Expression(in1);
            end
            if ~isnumeric(in2)
                in2 = BLOM_Expression(in2);
            end
            out = mrdivide(in1, in2);
        end
        
        function out = rdivide(in1, in2)
            if ~isnumeric(in1)
                in1 = BLOM_Expression(in1);
            end
            if ~isnumeric(in2)
                in2 = BLOM_Expression(in2);
            end
            out = rdivide(in1, in2);
        end
        
        function out = mpower(in1, in2)
            if ~isnumeric(in2)
                error('variables as exponents is not supported')
            elseif numel(in2) > 1
                error('non-scalar powers should use .^')
            elseif min(size(in1.idx)) > 1
                error('powers of matrix variables is not supported')
            elseif numel(in1.idx) > 1
                error('powers of vector variables should use .^')
            end
            out = BLOM_Expression(in1);
            out.Pt = in2*out.Pt;
            out.specialFunction = any(ismember(in2, ...
                BLOM_FunctionCode('all_codes')));
        end
        
        function out = power(in1, in2)
            if ~isnumeric(in2)
                error('variables as exponents is not supported')
            end
            size1 = size(in1.idx);
            size2 = size(in2);
            if max(size2) == 1
                % non-scalar variable to scalar power
                in2 = in2*ones(size1);
            elseif max(size1) > 1 && ~isequal(size1, size2)
                error('inputs must be scalars, or the same size as each other')
            end
            % scalar variable to non-scalar power, or non-scalar variable
            % to non-scalar power are both handled by sparse()
            out = BLOM_Expression(in1.problem, sparse(in1.idx(:), ...
                    1:numel(in2), in2(:), numel(in1.problem.x), ...
                    numel(in2)), speye(numel(in2)));
        end
        
        function out = ctranspose(in1)
            % transpose idx & value arrays
            out = in1;
            out.idx   = ctranspose(out.idx);
            out.value = ctranspose(out.value);
        end
        
        function out = transpose(in1)
            % transpose idx & value arrays
            out = in1;
            out.idx   = transpose(out.idx);
            out.value = transpose(out.value);
        end
        
        function out = dot(in1, in2)
            if ~isnumeric(in1)
                in1 = BLOM_Expression(in1);
            end
            if ~isnumeric(in2)
                in2 = BLOM_Expression(in2);
            end
            out = dot(in1, in2);
        end
        
        function out = prod(in1, dim)
            size1 = size(in1.idx);
            if nargin < 2 && min(size1) == 1
                % product of vector elements
                out = BLOM_Expression(in1.problem, sparse(in1.idx, 1, 1, ...
                    numel(in1.problem.x), 1), speye(1), false);
            elseif nargin < 2 || isequal(dim, 1)
                % product of matrix elements along columns
                
                
            elseif ~isnumeric(dim) || numel(dim) > 1
                error('dimension input must be a scalar constant')
            elseif isequal(dim, 2)
                % product of matrix elements along rows
                
            else
                error('multidimensional variables not supported')
            end
        end
        
        function out = sum(in1, dim)
            
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
        
        % norm
        % horzcat
        % vertcat
        % repmat
        % reshape?
        % subsasgn
        % le (output BLOM_Constraint object)
        % ge (output BLOM_Constraint object)
        % eq (output BLOM_Constraint object)
        
        % other trig, etc functions either derived or implicitly?
    end
end
