function data = BLOM_LoadDataFile(filename, formatstr)
% This function loads a data file into a Matlab variable, in the format
%   specified by the formatstr argument
% Input arguments:
% filename - the file name to save to
% formatstr - a string indicating the data format (case insensitive),
%   where the valid choices are:
%   tripletmat_binary for a sparse matrix in binary triplet [row col val] format
%   tripletmat_ascii  for a sparse matrix in ascii triplet [row col val] format
%   petscmat_binary   for a sparse matrix in petsc binary compressed sparse
%                     row format, can also use just petscmat
%   sparsevec_binary  for a sparse column vector in binary [row val] format
%   sparsevec_ascii   for a sparse column vector in ascii [row val] format
%   densevec_binary   for a dense vector in binary format
%   densevec_ascii    for a dense vector in ascii format
%   If formatstr is not provided, then the extension of the filename will
%   be used instead. If binary vs ascii is not specified, then the content
%   of the file will be examined and a heuristic will be used to guess.
% Output arguments:
% data - the loaded Matlab data

% NOTE: should probably rearrange tripletmat_binary and sparsevec_binary to
% group indices together (will need a header to indicate dimensions), since
% fwrite with a skip argument is slower than writing in contiguous chunks

if nargin == 1 || isempty(formatstr)
    % if formatstr not provided, then use extension of filename
    [dirname, basename, formatstr] = fileparts(filename);
    if ~isempty(formatstr) && formatstr(1) == '.'
        % remove initial .
        formatstr = formatstr(2:end);
    end
end
if strcmpi(formatstr, 'petscmat')
    % petscmat is always binary format
    formatstr = [formatstr '_binary'];
end

fid = fopen(filename, 'r');
if ~strcmpi(formatstr, 'txt') && ~strcmpi(formatstr(max(1,end-4):end), 'ascii') && ...
        ~strcmpi(formatstr(max(1,end-5):end), 'binary')
    % examine first 10000 bytes of file and use heuristic (like svn does)
    % to guess whether file is in ascii or binary format
    bytes = fread(fid, 10000, 'uint8=>uint8');
    if any(bytes == 0)
        % contains a byte equal to 0 (null character)
        formatstr = [formatstr '_binary'];
    else
        ascii_printable_chars = (bytes == 9) | (bytes == 10) | ...
            (bytes == 13) | (bytes > 31 & bytes < 127);
        if nnz(~ascii_printable_chars) > 0.15*length(bytes)
            % more than 15% of first 10000 bytes are not ascii printable characters
            formatstr = [formatstr '_binary'];
        else
            formatstr = [formatstr '_ascii'];
        end
    end
    frewind(fid)
end

switch lower(formatstr)
    case 'tripletmat_binary'
        % skip applies after reading each element for fread
        rows = fread(fid, inf, 'int', 12); % skip 12 bytes between row entries
        fseek(fid, 4, -1); % 4 bytes from start of file
        cols = fread(fid, inf, 'int', 12); % skip 12 bytes between col entries
        fseek(fid, 8, -1); % 8 bytes from start of file
        vals = fread(fid, inf, 'double', 8); % skip 8 bytes between val entries
    case 'tripletmat_ascii'
        tripletdata = fscanf(fid, '%g', [3 inf])';
        rows = tripletdata(:,1);
        cols = tripletdata(:,2);
        vals = tripletdata(:,3);
    case 'petscmat_binary'
        % first 4 entries of file are an integer flag MAT_FILE_CLASSID, the
        % number of rows, the number of columns, and the number of nonzeros
        header = fread(fid, 4, 'int');
        % next num_rows entries are the number of nonzeros per row
        nnz_per_row = fread(fid, header(2), 'int');
        nnz_prev_rows = [0; cumsum(nnz_per_row)]; % cumulative sum
        rows = zeros(header(4), 1); % preallocate
        for i=1:header(2)
            rows(nnz_prev_rows(i) + 1 : nnz_prev_rows(i + 1)) = i;
        end
        % next nnz entries are the 0-based column indices of the nonzeros
        cols = fread(fid, header(4), 'int') + 1;
        % last nnz entries are the nonzero values
        vals = fread(fid, header(4), 'double');
        if max(rows) ~= header(2) || max(cols) ~= header(3)
            % append an explicit zero to make sure matrix is created with
            % correct dimensions
            rows = [rows; header(2)];
            cols = [cols; header(3)];
            vals = [vals; 0];
        end
    case 'sparsevec_binary'
        % skip applies after reading each element for fread
        rows = fread(fid, inf, 'int', 8); % skip 8 bytes between ind entries
        fseek(fid, 4, -1); % 4 bytes from start of file
        vals = fread(fid, inf, 'double', 4); % skip 4 bytes between val entries
        cols = 1;
        if rows(1) == 0
            cols = [0; cols]; % insert cols(1) = 0 to remove below
        end
    case 'sparsevec_ascii'
        indsvals = load(filename);
        rows = indsvals(:,1);
        vals = indsvals(:,2);
        cols = 1;
        if rows(1) == 0
            cols = [0; cols]; % insert cols(1) = 0 to remove below
        end
    otherwise
        if strcmpi(formatstr, 'txt') || strcmpi(formatstr(max(1,end-4):end), 'ascii')
            try
                % try using load first, it can handle many cases
                data = load(filename);
            catch
                % if load doesn't work, just return the string contents
                data = fread(fid, [1 inf], 'char=>char');
            end
        else
            data = fread(fid, inf, 'double');
        end
end
fclose(fid);

if strncmpi(formatstr, 'tripletmat', 10) || strncmpi(formatstr, 'petscmat', 8) || ...
        strncmpi(formatstr, 'sparsevec', 9)
    if rows(1) == 0 && cols(1) == 0
        % 0, 0 position holds default_value if it is nonzero
        default_value = vals(1);
        rows = rows(2:end);
        cols = cols(2:end);
        vals = vals(2:end);
        data = default_value * ones(max(rows), max(cols)); % preallocate
        data(sub2ind(size(data), rows, cols)) = vals;
    else
        data = sparse(rows, cols, vals);
    end
end
