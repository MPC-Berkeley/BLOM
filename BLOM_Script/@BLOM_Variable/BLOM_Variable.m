% The BLOM_Variable class represents optimization variables according to
% their parent BLOM_Problem object and the index within that problem.
% It should NOT be used directly to create new variables! 
% Use the newVariable method of the BLOM_Problem class to do that.

classdef BLOM_Variable < handle
    properties
        problem % handle of parent problem
        idx % index within parent problem
        value % numerical value (set to nan if no values known)
    end
    
    methods
        function var = BLOM_Variable(problem, idx, value)
            var.problem = problem;
            var.idx = idx;
            if nargin < 3
                var.value = nan(size(idx));
            else
                if ~isequal(size(idx), size(value))
                    error('idx and value must be same size')
                end
                var.value = value;
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
                    sparse(numel(in2.problem.lb), 1), in1(:), false);
            elseif isnumeric(in2)
                in2 = BLOM_Expression(in1.problem, ...
                    sparse(numel(in1.problem.lb), 1), in2(:), false);
            end
            out = plus(BLOM_Expression(in1), BLOM_Expression(in2));
        end
        
        function out = minus(in1, in2)
            if isnumeric(in1)
                in1 = BLOM_Expression(in2.problem, ...
                    sparse(numel(in2.problem.lb), 1), in1(:), false);
            elseif isnumeric(in2)
                in2 = BLOM_Expression(in1.problem, ...
                    sparse(numel(in1.problem.lb), 1), in2(:), false);
            end
            out = minus(BLOM_Expression(in1), BLOM_Expression(in2));
        end
        
        function out = uminus(in1) % -x
            out = uminus(BLOM_Expression(in1));
        end
        
        function out = mtimes(in1, in2)
            if isnumeric(in1)
                out = mtimes(in1, BLOM_Expression(in2));
            elseif isnumeric(in2)
                out = mtimes(BLOM_Expression(in1), in2);
            else
                out = mtimes(BLOM_Expression(in1), BLOM_Expression(in2));
            end
        end
        
        function out = times(in1, in2)
            if isnumeric(in1)
                out = times(in1, BLOM_Expression(in2));
            elseif isnumeric(in2)
                out = times(BLOM_Expression(in1), in2);
            else
                out = times(BLOM_Expression(in1), BLOM_Expression(in2));
            end
        end
        
        function out = mrdivide(in1, in2)
            if isnumeric(in1)
                out = mrdivide(in1, BLOM_Expression(in2));
            elseif isnumeric(in2)
                out = mrdivide(BLOM_Expression(in1), in2);
            else
                out = mrdivide(BLOM_Expression(in1), BLOM_Expression(in2));
            end
        end
        
        function out = rdivide(in1, in2)
            if isnumeric(in1)
                out = rdivide(in1, BLOM_Expression(in2));
            elseif isnumeric(in2)
                out = rdivide(BLOM_Expression(in1), in2);
            else
                out = rdivide(BLOM_Expression(in1), BLOM_Expression(in2));
            end
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
                    1:numel(in2), in2(:), numel(in1.problem.lb), ...
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
            if isnumeric(in1)
                size1 = size(in1);
                size2 = size(in2.idx);
                if min(size1) > 1 || min(size2) > 1
                    error('both inputs must be vectors')
                elseif max(size1) ~= max(size2)
                    error('both inputs must be the same length')
                end
                % FINISH HERE
            elseif isnumeric(in2)
                size1 = size(in1.idx);
                size2 = size(in2);
                if min(size1) > 1 || min(size2) > 1
                    error('both inputs must be vectors')
                elseif max(size1) ~= max(size2)
                    error('both inputs must be the same length')
                end
                % FINISH HERE
            else
                out = dot(BLOM_Expression(in1), BLOM_Expression(in2));
            end
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
