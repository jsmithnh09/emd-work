function y = db2mag(x)
% MAG2DB performs decibel to magnitude conversion.
%
%   Y = DB2MAG(X)
%
% See Also: MAG2DB

y = 10.^(x ./ 20);

end