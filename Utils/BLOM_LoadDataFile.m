function data = BLOM_LoadDataFile(filename, formatstr)
% This function loads a data file into a Matlab variable, in the format
%   specified by the formatstr argument
% Input arguments:
% filename - the file name to save to
% formatstr - a string indicating the data format (case insensitive),
%   where the valid choices are:
%   tripletmat_binary for a sparse matrix in binary triplet [row col val] format
%   tripletmat_ascii  for a sparse matrix in ascii triplet [row col val] format
%   sparsevec_binary  for a sparse column vector in binary [row val] format
%   sparsevec_ascii   for a sparse column vector in ascii [row val] format
%   densevec_binary   for a dense vector in binary format
%   densevec_ascii    for a dense vector in ascii format
%   If formatstr is not provided, then the extension of the filename will
%   be used instead. If binary vs ascii is not specified, then the content
%   of the file will be examined and a heuristic will be used to guess.
% Output arguments:
% data - the loaded Matlab data

if nargin == 1 || isempty(formatstr)
    % if formatstr not provided, then use extension of filename
    [dirname, basename, formatstr] = fileparts(filename);
    if ~isempty(formatstr) && formatstr(1) == '.'
        % remove initial .
        formatstr = formatstr(2:end);
    end
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
        ascii_printable_chars = [9, 10, 13, 32:126];
        num_unprintable = nnz(~any(repmat(bytes, 1, length(ascii_printable_chars)) == ...
            repmat(ascii_printable_chars, length(bytes), 1), 2));
        if num_unprintable > 0.15*length(bytes)
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

if strncmpi(formatstr, 'tripletmat', 10) || strncmpi(formatstr, 'sparsevec', 9)
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
