function s = nan2hex(x)
% Very similar function to Matlab's num2hex, but accurately outputs the bit
% contents of NaN values. Some inconsistency may occur for signaling NaNs.
% This function is not properly vectorized yet, but that can be fixed.

x = typecast(x,'uint64');
if ~bitget(x,52) && all(bitget(x,53:63)) && any(bitget(x,1:51))
   warning(['Signaling NaN detected, leading mantissa bit may be modified. ' ...
      'Quiet NaNs are in the ranges 7ff8000000000000 to 7fffffffffffffff ' ...
      'or fff8000000000000 to ffffffffffffffff.'])
end
list1 = uint64(2.^(0:4:60));
d = '0123456789abcdef';
hexdigits = fliplr(diff([mod(x,list1) x])./list1);
s = d(hexdigits+1);
