function [lower, upper] = gcmci(P, N, Z)
% GGD GCM confidence interval from GCM search method.
%
%   [LOWER, UPPER] = GGDCI(P, N, Z)
%       P (scalar) is the estimate
%       N (scalar) are the number of samples
%       Z (scalar) is the z-score (default: 2.576)
%
% See Also: GCM_SEARCH

if (nargin < 3)
    Z = 2.576;
end
lower = (coth(acoth(sqrt(P + 1)) + ((1/(sqrt(2*N)))*Z))^2) - 1;
upper = (coth(acoth(sqrt(P + 1)) - ((1/(sqrt(2*N)))*Z))^2) - 1;

end