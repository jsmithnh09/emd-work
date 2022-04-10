function y = ggdrnd(mu, sigma, beta, dims)
%GGDRND generates random samples from the GGD.
%   Y = RANDGGD(MU, SIGMA, BETA, DIMS) returns an array of random numbers
%   from the GGD distribution with location parameter MU, variation 
%   parameter SIGMA, and shape parameter BETA.
%
%   The scale parameter alpha is computed based on the variance equation,
%   since variance and shape are provided.
%
%   Y = GGDRND(MU, SIGMA, BETA, N) returns an N-by-1 array.
%   Y = GGDRND(MU, SIGMA, BETA, [M, N, ... P]) returns an M-by-N-etc.
%       matrix.
% 
% See Also GAMRND, GAMMA

% References:
%   [1] Gonzalez-Farias, G., Molina, J. A. D., & Rodríguez-Dagnino, 
%       R. M. (2009). Efficiency of the approximated shape parameter 
%       estimator in the generalized Gaussian distribution. IEEE 
%       Transactions on Vehicular Technology, 58(8), 4214-4223.

if (isinf(beta) || (beta <= 0))
    error('ggdrnd:InvalidShape', ...
        'Shape parameter outside (0, Inf) range.');
elseif (numel(dims) == 1)
    dims = [dims(1), 1];
end

alpha = sigma * sqrt(gamma(1/beta) / gamma(3/beta)); % Scale
b = 2*(rand(dims) < 0.5) - 1; % Bernoulli
g = gamrnd(1/beta, 1, dims).^(1/beta);
y = mu + (1/sqrt(alpha)) .* g .* b;
y = y(:);

end
