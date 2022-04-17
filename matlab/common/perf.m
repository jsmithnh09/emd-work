function [snr, mse, rmse] = perf(X, Y)
% PERF determines the filter performance.
%
%   [SNR, MSE, RMSE] = PERFORMANCE(X, Y)
%       X (vector) is the input signal.
%       Y (vector) is the output signal.
%
% See Also: IMMSE

X = X(:);
Y = Y(:);
if (length(X) ~= length(Y))
    error('perf::mismatched vector lengths.');
end
snr = 10*log10(sum(Y.^2) / sum((X-Y).^2));
mse = (norm(X-Y, 2).^2) / numel(X);
rmse = sqrt(mse);

end % performance
