function out = power(in1, in2)

if ~isnumeric(in2)
    error('variables as exponents is not supported')
end
size1 = [size(in1.K, 1), 1];
size2 = size(in2);
if size2(1) == 1 && size2(2) > 1
    warning(['all BLOM_Expression objects are considered vectors, ' ...
        'converting second input of power to column vector'])
    in2 = in2(:);
    size2 = size(in2);
end
if max(size1) > 1 && max(size2) > 1 && ~isequal(size1, size2)
    error('inputs must be scalars, or the same size as each other')
end
specialFunction = any(ismember(in2, BLOM_FunctionCode('all_codes')));
in1 = in1.removeUnusedTerms;
in1_Kt = in1.K'; % save transpose because getting columns of a sparse
% matrix is faster than getting rows
terms_per_row = full(sum(spones(in1_Kt), 1))';
if ~specialFunction && ~in1.specialFunction && all(terms_per_row <= 1)
    % each output row is a single-term expression with no special functions,
    % so can just multiply all exponents in Pt and take powers of K's nonzeros
    if max(size2) == 1
        % non-scalar expression to scalar power
        if in2 == 0
            % output is a vector of 1's
            out = BLOM_Expression(in1.problem, ...
                sparse(numel(in1.problem.lb), 1), ones(size1), false);
        else
            out = in1;
            out.Pt = in2*out.Pt;
            out.K(out.K ~= 0) = (out.K(out.K ~= 0)).^in2;
        end
    elseif max(size1) == 1
        % scalar expression to non-scalar power
        newPt = kron(in2', in1.Pt);
        newK = sparse(1:numel(in2), 1:numel(in2), (in1.K).^in2);
        out = BLOM_Expression(in1.problem, newPt, newK, false);
        out.auxPt = in1.auxPt;
        out.auxK  = in1.auxK;
    else
        % non-scalar expression to non-scalar power
        newPt = cell(size2);
        newK_diags = zeros(size2);
        for i = 1:numel(in2) % should vectorize this loop
            if in2(i) == 0
                % this output is a constant 1
                newPt{i} = sparse(numel(in1.problem.lb), 1);
                newK_diags(i) = 1;
            elseif terms_per_row(i) == 0
                % this output is a constant 0
                newPt{i} = sparse(numel(in1.problem.lb), 1);
                newK_diags(i) = 0;
            else
                K_row = in1_Kt(:, i);
                newPt{i} = in1.Pt(:, K_row ~= 0) * in2(i);
                newK_diags(i) = nonzeros(K_row)^in2(i);
            end
        end
        newK = sparse(1:numel(in2), 1:numel(in2), newK_diags);
        out = BLOM_Expression(in1.problem, [newPt{:}], newK, false);
        out.auxPt = in1.auxPt;
        out.auxK  = in1.auxK;
    end
else
    % have to introduce new variable to take power of a multi-term
    % expression or special function, or special function of any expression
    newvar = in1.problem.newVariable(size1);
    out = power(newvar, in2);
    % set aux Pt and K in output expression for newvar == in1 constraint
    aux = in1 - newvar;
    out.auxPt = [aux.Pt, aux.auxPt];
    out.auxK = blkdiag(aux.K, aux.auxK);
end
