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
%   sparsevec_binary  for a sparse vector in binary [ind val] format
%   sparsevec_ascii   for a sparse vector in ascii [ind val] format
%   densevec_binary   for a dense vector in binary format
%   densevec_ascii    for a dense vector in ascii format
% default_value - value used for unspecified elements in sparse formats,
%   assumed to be 0 if not given (use inf for ub, -inf for lb)

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
        [rows, cols, vals] = find(data);
    else
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
    if ~any(nondefault_elements(:,end)) || ~any(nondefault_elements(end,:))
        % append a default_value if last row or last column is empty
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
