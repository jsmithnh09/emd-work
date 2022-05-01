function y = db2pow(x)
% MAG2DB performs decibel to power conversion.
%
%   Y = DB2POW(X)
%
% See Also: DB2MAG

y = 10.^(x ./ 10);

end
