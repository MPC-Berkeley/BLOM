function x = hex2nan(s)
% Very similar function to Matlab's hex2num, but does not change the bit
% contents of NaN values. Some inconsistency may occur for signaling NaNs.
% This function is not properly vectorized yet, but that can be fixed.

[row,col] = size(s);
blanks = find(s==' '); % Find the blanks at the end
if ~isempty(blanks)
    s(blanks) = '0';
end % Zero pad the shorter hex numbers.

% Convert characters to numeric digits.
% More than 16 characters are truncated.
d = zeros(row,16);
d(:,1:col) = abs(lower(s)) - '0';
d = d + ('0'+10-'a').*(d>9);

x = sum(uint64(d).*uint64(2.^(60:-4:0)),'native');
if ~bitget(x,52) && all(bitget(x,53:63)) && any(bitget(x,1:51))
   warning(['Signaling NaN detected, leading mantissa bit may be modified. ' ...
      'Quiet NaNs are in the ranges 7ff8000000000000 to 7fffffffffffffff ' ...
      'or fff8000000000000 to ffffffffffffffff.'])
end
x = typecast(x,'double');
