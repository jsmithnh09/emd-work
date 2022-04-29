function kcoeff = ggdcsk(X)
% GGDCSK GGD kurtosis coefficient.
%
%   KCOEFF = GGDCSK(X)
%       X (double) is a random variable.
%       KCOEFF (scalar) is the kurtosis of the sample.
%
% See Also: GGDKURT

X = X - mean(X);
N = length(X);
kcoeff = N * (sum(X.^4)/sum(X.^2)^2) - 3;

end