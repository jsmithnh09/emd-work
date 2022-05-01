function y = mag2db(x)
% MAG2DB performs magnitude to decibel conversion.
%
%   Y = MAG2DB(X)
%
% See Also: DB2MAG

y = 20.*log10(x);

end