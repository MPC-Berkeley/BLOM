function BLOM_SaveDataFile(data, filename, formatstr, default_value)
% This function saves a Matlab variable to a data file, in the format
%   specified by the formatstr argument
% Input arguments:
% data - the Matlab data to save
% filename - the file name to save to
% formatstr - a string indicating the data format (case insensitive),
%   where the valid choices are:
%   tripletmat_binary for a sparse matrix in binary triplet [row col val] format
%   tripletmat_ascii  for a sparse matrix in ascii triplet [row col val] format
%   petscmat_binary   for a sparse matrix in petsc binary compressed sparse
%                     row format, can also use just petscmat
%   sparsevec_binary  for a sparse vector in binary [ind val] format
%   sparsevec_ascii   for a sparse vector in ascii [ind val] format
%   densevec_binary   for a dense vector in binary format
%   densevec_ascii    for a dense vector in ascii format
% default_value - value used for unspecified elements in sparse formats,
%   assumed to be 0 if not given (use inf for ub, -inf for lb)

% NOTE: should probably rearrange tripletmat_binary and sparsevec_binary to
% group indices together (will need a header to indicate dimensions), since
% fwrite with a skip argument is slower than writing in contiguous chunks

if ~strncmpi(formatstr, 'densevec', 8)
    if nargin < 4
        default_value = 0;
    end
    if isnan(default_value)
        nondefault_elements = ~isnan(data);
    else
        nondefault_elements = (data ~= default_value);
    end
    % find non-default rows, cols, vals for sparse mat/vec cases
    if default_value == 0
        if strncmpi(formatstr, 'petscmat', 8)
            % row-major order for petsc format
            [cols, rows, vals] = find(data');
        else
            % default column-major order for others
            [rows, cols, vals] = find(data);
        end
    else
        if strncmpi(formatstr, 'petscmat', 8)
            error('petscmat format must use default_value=0')
        end
        [rows, cols] = find(nondefault_elements);
        vals = full(data(sub2ind(size(data), rows, cols)));
    end
    % force rows, cols, vals to be column vectors in case data is a row vector
    rows = reshape(rows, [], 1);
    cols = reshape(cols, [], 1);
    vals = reshape(vals, [], 1);
    if default_value ~= 0
        % save default_value in the 0, 0 position if it is nonzero
        rows = [0; rows];
        cols = [0; cols];
        vals = [default_value; vals];
    end
    if (~any(nondefault_elements(:,end)) || ~any(nondefault_elements(end,:))) ...
            && ~strncmpi(formatstr, 'petscmat', 8)
        % append a default_value if last row or last column is empty
        % for triplet or sparse vector formats
        rows = [rows; size(data, 1)];
        cols = [cols; size(data, 2)];
        vals = [vals; default_value];
    end
    if strncmpi(formatstr, 'sparsevec', 9)
        if size(data,2) == 1
            % data is a column vector
            inds = rows;
        else
            % data is a row vector
            inds = cols;
        end
    end
end

fid = fopen(filename, 'w');
switch lower(formatstr)
    case 'tripletmat_binary'
        % skip applies before writing each element for fwrite,
        % so need to write first element, then rest with skipping
        fwrite(fid, [rows(1); cols(1)], 'int');
        fwrite(fid, [rows(2:end), cols(2:end)]', '2*int', 8); % skip 8 bytes between pairs of index entries
        frewind(fid);
        fwrite(fid, vals, 'double', 8); % skip 8 bytes between val entries
    case 'tripletmat_ascii'
        fprintf(fid, '%d %d %.17g\n', [rows, cols, vals]');
    case {'petscmat_binary', 'petscmat'}
        % compressed sparse row format with 0-based indices
        [num_rows, num_cols] = size(data);
        if min(size(data)) == 1
            nnz_per_row = full(data' ~= 0);
        else
            nnz_per_row = full(sum(data' ~= 0, 1));
        end
        % first 4 entries of file are an integer flag MAT_FILE_CLASSID, the
        % number of rows, the number of columns, and the number of nonzeros
        fwrite(fid, [1211216, num_rows, num_cols, nnz(data)], 'int');
        % next num_rows entries are the number of nonzeros per row
        fwrite(fid, nnz_per_row, 'int');
        % next nnz(data) entries are the 0-based column indices of the nonzeros
        fwrite(fid, cols - 1, 'int');
        % last nnz(data) entries are the nonzero values
        fwrite(fid, vals, 'double');
    case 'sparsevec_binary'
        % skip applies before writing each element for fwrite,
        % so need to write first element, then rest with skipping
        fwrite(fid, inds(1), 'int');
        fwrite(fid, inds(2:end), 'int', 8); % skip 8 bytes between ind entries
        frewind(fid);
        fwrite(fid, vals, 'double', 4); % skip 4 bytes between val entries
    case 'sparsevec_ascii'
        fprintf(fid, '%d %.17g\n', [inds, vals]');
    case 'densevec_binary'
        fwrite(fid, data, 'double');
    case 'densevec_ascii'
        fprintf(fid, '%.17g\n', data);
    otherwise
        fclose(fid);
        error(['Unrecognized formatstr ' formatstr])
end
fclose(fid);
