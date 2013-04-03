function out = times(in1, in2)

if isnumeric(in1)
    size1 = size(in1);
    if size1(2) > 1
        warning(['all BLOM_Expression objects are considered vectors, ' ...
            'converting first input of times to column vector'])
        in1 = in1(:);
        size1 = size(in1);
    end
    problem = in2.problem;
else
    % convert to BLOM_Expression (refreshes size of in1.Pt)
    in1 = BLOM_Expression(in1);
    size1 = [size(in1.K, 1), 1];
    problem = in1.problem;
end
if isnumeric(in2)
    size2 = size(in2);
    if size2(2) > 1
        warning(['all BLOM_Expression objects are considered vectors, ' ...
            'converting second input of times to column vector'])
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
if ~isequal(size1, size2) && max(size1) > 1 && max(size2) > 1
    error('dimensions not compatible for elementwise multiplication')
end
if isnumeric(in1)
    if max(size2) == 1
        % constant vector times scalar expression
        out = in2;
        out.K = kron(in1, out.K);
    else
        % constant vector times vector expression
        % (or constant scalar times vector expression)
        out = in2;
        out.K = sparse(1:numel(in1), 1:numel(in1), in1) * out.K;
    end
elseif isnumeric(in2)
    if max(size1) == 1
        % scalar expression times constant vector
        out = in1;
        out.K = kron(in2, out.K);
    else
        % vector expression times constant vector
        % (or vector expression times constant scalar)
        out = in1;
        out.K = sparse(1:numel(in2), 1:numel(in2), in2) * out.K;
    end
else % expression times expression
    in1 = in1.removeUnusedTerms;
    in2 = in2.removeUnusedTerms;
    in1_Kt = in1.K';
    in1_terms_per_row = full(sum(spones(in1_Kt), 1))';
    in2_Kt = in2.K';
    in2_terms_per_row = full(sum(spones(in2_Kt), 1))';
    newvar = [];
    if in1.specialFunction || in2.specialFunction || ...
            (any(in1_terms_per_row > 1) && any(in2_terms_per_row > 1))
        % if either expression has special functions or both expressions have
        % multiple terms, introduce a variable for the scalar input (if any)
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
        elseif in2.specialFunction && ~in1.specialFunction
            % both inputs are vectors, but only in2 has special functions
            % so introduce vector variable for in2
            newvar = problem.newVariable(size2);
            newPt = cell(size1);
            newK = cell(size1);
            for i = 1:max(size1) % should vectorize this loop
                terms1 = (in1_Kt(:, i) ~= 0);
                if ~any(terms1)
                    % this element of in1 is zero
                    newPt{i} = [in1.Pt(:, terms1); sparse(max(size2), 0)];
                    newK{i} = sparse(1, 0);
                else
                    % multiply every term of this element by newvar(i)
                    newPt{i} = [in1.Pt(:, terms1); sparse(i, ...
                        1:nnz(terms1), 1, max(size2), nnz(terms1))];
                    newK{i} = in1_Kt(terms1, i)';
                end
            end
            out = BLOM_Expression(problem, [newPt{:}], blkdiag(newK{:}), false);
            % set aux Pt and K in output expression for newvar == in2 constraint
            aux = in2 - newvar;
        else
            % both inputs are vectors, introduce vector variable for in1
            newvar = problem.newVariable(size1);
            newPt = cell(size2);
            newK = cell(size2);
            for i = 1:max(size2) % should vectorize this loop
                terms2 = (in2_Kt(:, i) ~= 0);
                if ~any(terms2)
                    % this element of in2 is zero
                    newPt{i} = [in2.Pt(:, terms2); sparse(max(size1), 0)];
                    newK{i} = sparse(1, 0);
                else
                    % multiply every term of this element by newvar(i)
                    newPt{i} = [in2.Pt(:, terms2); sparse(i, ...
                        1:nnz(terms2), 1, max(size1), nnz(terms2))];
                    newK{i} = in2_Kt(terms2, i)';
                end
            end
            out = BLOM_Expression(problem, [newPt{:}], blkdiag(newK{:}));
            % set aux Pt and K in output expression for newvar == in1 constraint
            aux = in1 - newvar;
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
    elseif all(in1_terms_per_row <= 1)
        % in1 is single term per vector element, in2 is vector multi-term
        newPt = cell(size1);
        newK = cell(size1);
        for i = 1:max(size1) % should vectorize this loop
            term1 = (in1_Kt(:, i) ~= 0);
            terms2 = (in2_Kt(:, i) ~= 0);
            if ~any(term1) || ~any(terms2)
                % this element of in1 or in2 is zero
                newPt{i} = in1.Pt(:, term1);
                newK{i} = sparse(1, 0);
            else
                % multiply all the terms that are used in element i of in2
                % by the single term in element i of in1
                newPt{i} = in2.Pt(:, terms2) + repmat(in1.Pt(:, term1), ...
                    1, nnz(terms2));
                newK{i} = in2_Kt(terms2, i)' * in1_Kt(term1, i);
            end
        end
        out = BLOM_Expression(problem, [newPt{:}], blkdiag(newK{:}), false);
    else
        % in2 is single term per vector element, in1 is vector multi-term
        newPt = cell(size1);
        newK = cell(size1);
        for i = 1:max(size1) % should vectorize this loop
            terms1 = (in1_Kt(:, i) ~= 0);
            term2 = (in2_Kt(:, i) ~= 0);
            if ~any(terms1) || ~any(term2)
                % this element of in1 or in2 is zero
                newPt{i} = in1.Pt(:, terms1);
                newK{i} = sparse(1, 0);
            else
                % multiply all the terms that are used in element i of in1
                % by the single term in element i of in2
                newPt{i} = in1.Pt(:, terms1) + repmat(in2.Pt(:, term2), ...
                    1, nnz(terms1));
                newK{i} = in1_Kt(terms1, i)' * in2_Kt(term2, i);
            end
        end
        out = BLOM_Expression(problem, [newPt{:}], blkdiag(newK{:}), false);
    end
    out.auxPt = [in1.auxPt, in2.auxPt];
    out.auxK = blkdiag(in1.auxK, in2.auxK);
    if ~isempty(newvar) % introduced new variable, add to auxPt and auxK
        out.auxPt = [blkdiag(out.auxPt, ...
            sparse(max(size(newvar.idx)), 0)), aux.Pt];
        out.auxK = blkdiag(out.auxK, aux.K);
    end
end
