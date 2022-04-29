function kcoeff = ggdkurt(B)
% GGDKURT GGD theoretical kurtosis.
%
%   KCOEFF = GGDKURT(X)
%       X (double) is a random variable.
%       KCOEFF (scalar) is the kurtosis of the sample.
%
% See Also: GGDKURT

kcoeff = ((gamma(5./B) .* gamma(1./B)) ./ (gamma(3./B).^2)) - 3;

end