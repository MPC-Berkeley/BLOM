function out = plus(in1, in2)

if isnumeric(in1)
    size1 = size(in1);
    if size1(2) > 1
        size1 = [prod(size1), 1];
        warning(['all BLOM_Expression objects are considered vectors, ' ...
                'converting first input of plus to column vector'])
    end
    in1 = BLOM_Expression(in2.problem, ...
        sparse(numel(in2.problem.x), 1), in1(:), false);
else
    % convert to BLOM_Expression (refreshes size of in1.Pt)
    in1 = BLOM_Expression(in1);
    size1 = [size(in1.K, 1), 1];
end
if isnumeric(in2)
    size2 = size(in2);
    if size2(2) > 1
        size2 = [prod(size2), 1];
        warning(['all BLOM_Expression objects are considered vectors, ' ...
                'converting second input of plus to column vector'])
    end
    in2 = BLOM_Expression(in1.problem, ...
        sparse(numel(in1.problem.x), 1), in2(:), false);
else
    % convert to BLOM_Expression (refreshes size of in2.Pt)
    in2 = BLOM_Expression(in2);
    size2 = [size(in2.K, 1), 1];
end
if ~isequal(in1.problem, in2.problem)
    error('cannot combine expressions from different problems')
elseif max(size1) > 1 && max(size2) > 1 && ~isequal(size1, size2)
    error('inputs must be scalars, or the same size as each other')
end
% concatenate K's, repeated if one input is a scalar and the other isn't
if max(size1) > max(size2)
    newK = [in1.K, repmat(in2.K, max(size1), 1)];
elseif max(size2) > max(size1)
    newK = [repmat(in1.K, max(size2), 1), in2.K];
else
    newK = [in1.K, in2.K];
end
out = BLOM_Expression(in1.problem, [in1.Pt, in2.Pt], newK, ...
    in1.specialFunction || in2.specialFunction);
out.auxPt = [in1.auxPt, in2.auxPt];
out.auxK = blkdiag(in1.auxK, in2.auxK);
