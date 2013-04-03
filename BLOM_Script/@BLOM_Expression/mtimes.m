function out = mtimes(in1, in2)

if isnumeric(in1)
    size1 = size(in1);
    problem = in2.problem;
else
    % convert to BLOM_Expression (refreshes size of in1.Pt)
    in1 = BLOM_Expression(in1);
    size1 = [size(in1.K, 1), 1];
    problem = in1.problem;
end
if isnumeric(in2)
    size2 = size(in2);
    if size2(2) > 1 && max(size1) > 1
        error(['all BLOM_Expression objects are considered vectors, ' ...
            'so matrix multiplication is only supported from the left'])
    elseif size2(2) > 1
        warning(['all BLOM_Expression objects are considered vectors, ' ...
            'converting second input of mtimes to column vector'])
        in2 = in2(:);
        size2 = size(in2);
    end
else
    % convert to BLOM_Expression (refreshes size of in2.Pt)
    in2 = BLOM_Expression(in2);
    size2 = [size(in2.K, 1), 1];
    if ~isequal(problem, in2.problem)
        error('cannot combine expressions from different problems')
    end
end
if size1(2) ~= size2(1) && max(size1) > 1 && max(size2) > 1
    error(['dimensions not compatible for matrix multiplication, ' ...
        'use .* for elementwise multiplication'])
end
if isnumeric(in1)
    if max(size2) == 1
        % constant matrix times scalar expression
        if size1(2) > 1
            warning(['all BLOM_Expression objects are considered vectors, ' ...
                'converting first input of mtimes to column vector'])
            in1 = in1(:);
            size1 = size(in1);
        end
        out = in2;
        out.K = kron(in1, out.K);
    else
        % constant matrix times vector expression
        % (or constant scalar times any expression)
        out = in2;
        out.K = in1*out.K;
    end
elseif isnumeric(in2)
    if max(size2) == 1
        % any expression times constant scalar
        out = in1;
        out.K = in2*out.K;
    elseif max(size1) == 1
        % scalar expression times constant vector
        out = in1;
        out.K = kron(in2, out.K);
    else
        error('this case should not be possible')
    end
else % expression times expression, and we know one (or both) of the
    % expressions must be a scalar for mtimes not to error
    in1 = in1.removeUnusedTerms;
    in2 = in2.removeUnusedTerms;
    in1_Kt = in1.K';
    in1_terms_per_row = full(sum(spones(in1_Kt), 1))';
    in2_Kt = in2.K';
    in2_terms_per_row = full(sum(spones(in2_Kt), 1))';
    newvar = [];
    if in1.specialFunction || in2.specialFunction || ...
            (any(in1_terms_per_row > 1) && any(in2_terms_per_row > 1))
        % if either expression has special functions or both expressions
        % have multiple terms, introduce a variable for the scalar input
        % can maybe be more careful about special functions, only
        % introducing variables if the same variable used in special
        % function is also present in other input
        if max(size1) == 1 && max(size2) == 1 && in2.specialFunction && ...
                ~in1.specialFunction
            % if both inputs are scalar and only in2 is a special function,
            % introduce variable for in2
            newvar = problem.newVariable(size2);
            out = in1;
            % multiply all terms by newvar
            out.Pt = [out.Pt; ones(1, size(out.Pt, 2))];
            % set aux Pt and K in output expression for newvar == in2 constraint
            aux = in2 - newvar;
        elseif max(size1) == 1
            newvar = problem.newVariable(size1);
            out = in2;
            % multiply all terms by newvar
            out.Pt = [out.Pt; ones(1, size(out.Pt, 2))];
            % set aux Pt and K in output expression for newvar == in1 constraint
            aux = in1 - newvar;
        elseif max(size2) == 1
            newvar = problem.newVariable(size2);
            out = in1;
            % multiply all terms by newvar
            out.Pt = [out.Pt; ones(1, size(out.Pt, 2))];
            % set aux Pt and K in output expression for newvar == in2 constraint
            aux = in2 - newvar;
        else
            error('this case should not be possible')
        end
    elseif all(in1_terms_per_row <= 1) && (max(size1) == 1)
        % in1 is single term and scalar
        out = in2;
        out.K = out.K * in1.K;
        out.Pt = out.Pt + repmat(in1.Pt, 1, size(out.Pt, 2));
    elseif all(in2_terms_per_row <= 1) && (max(size2) == 1)
        % in2 is single term and scalar
        out = in1;
        out.K = out.K * in2.K;
        out.Pt = out.Pt + repmat(in2.Pt, 1, size(out.Pt, 2));
    elseif all(in1_terms_per_row <= 1) && (max(size2) == 1)
        % in1 is single term per vector element, in2 is scalar multi-term
        newPt = cell(size1);
        newK = cell(size1);
        for i = 1:max(size1) % should vectorize this loop
            term1 = (in1_Kt(:, i) ~= 0);
            if ~any(term1)
                % this element of in1 is zero
                newPt{i} = in1.Pt(:, term1);
                newK{i} = sparse(1, 0);
            else
                newPt{i} = in2.Pt + repmat(in1.Pt(:, term1), 1, size(in2.Pt, 2));
                newK{i} = in2.K * in1_Kt(term1, i);
            end
        end
        out = BLOM_Expression(problem, [newPt{:}], blkdiag(newK{:}), false);
    elseif all(in2_terms_per_row <= 1) && (max(size1) == 1)
        % in2 is single term per vector element, in1 is scalar multi-term
        newPt = cell(size2);
        newK = cell(size2);
        for i = 1:max(size2) % should vectorize this loop
            term2 = (in2_Kt(:, i) ~= 0);
            if ~any(term2)
                % this element of in2 is zero
                newPt{i} = in2.Pt(:, term2);
                newK{i} = sparse(1, 0);
            else
                newPt{i} = in1.Pt + repmat(in2.Pt(:, term2), 1, size(in1.Pt, 2));
                newK{i} = in1.K * in2_Kt(term2, i);
            end
        end
        out = BLOM_Expression(problem, [newPt{:}], blkdiag(newK{:}), false);
    else
        error('this case should not be possible')
    end
    out.auxPt = [in1.auxPt, in2.auxPt];
    out.auxK = blkdiag(in1.auxK, in2.auxK);
    if ~isempty(newvar) % introduced new variable, add to auxPt and auxK
        out.auxPt = [blkdiag(out.auxPt, ...
            sparse(max(size(newvar.idx)), 0)), aux.Pt];
        out.auxK = blkdiag(out.auxK, aux.K);
    end
end
