function [estimate, total] = gcm_search(X, initial, trials)
% searches GCM function with Newton-Raphson iterative search [1].
%
%   [ESTIMATE, TOTAL] = GCM_SEARCH(X[, INITIAL, TRIALS])
%       X (vector) is a random sample.
%       INITIAL (scalar) is the first guess of the distribution shape.
%       (Default: m1/sqrt(m2), "m" denotes statistical moments)
%       TRIALS (scalar) are the number of iterations to go through.
%       (Default: 100)
%
% Input:
%		X		-	input data sequence.
%		initial - 	first guess of the true shape.
%		trials	-	number of trials.
%
% If INITIAL is left empty, the first iteration is estimated as
% m1/sqrt(m2), where m denotes statistical moments.
%
% [1] Song, Kai-Sheng. "A globally convergent and consistent method
%   for estimating the shape parameter of a generalized Gaussian distribution."
%   IEEE Transactions on Information Theory 52.2 (2006): 510-527.

% initial values.

N = length(X);
if (nargin < 3)
    trials = 100;
end
if (nargin < 2) || (isempty(initial))
    initial = ((1/N) * sum(abs(X))) / std(X);
end
[Zn, ZnP] = deal(zeros(1, trials-1));
a = [initial, zeros(1, trials-1)];

Zn(1) = ((1/N)*sum(abs(X).^(2*a(1)))) / (((1/N) .* sum(abs(X).^(a(1))) )^2) - (1 + a(1) );
ZnP(1) = ((((2/N).*sum(abs(X).^(2.*a(1)).*log(abs(X)))) .* (((1/N).*sum(abs(X).^(a(1)))).^2) ) / (((1/N).*sum(abs(X).^(a(1)))).^4)) - ...
    ((((1/N).*sum(abs(X).^(a(1)).*log(abs(X)))) .* ((1/N).*sum(abs(X).^(2*a(1)))) .* ((2/N).*sum(abs(X).^(a(1))))) / (((1/N).*sum(abs(X).^(a(1)))).^4) ) - 1;

for i = 2:trials
    % Z(theta).
    Zn(i-1) = ((1/N).*sum(abs(X).^(2*a(i-1)))) / (((1/N) .* sum(abs(X).^(a(i-1))) )^2) - (1 + a(i-1) );

    % Z(theta) prime derivative.
    ZnP(i-1) = ((((2/N).*sum(abs(X).^(2.*a(i-1)).*log(abs(X)))) .* (((1/N).*sum(abs(X).^(a(i-1)))).^2) ) / (((1/N).*sum(abs(X).^(a(i-1)))).^4)) - ...
        ((((1/N).*sum(abs(X).^(a(i-1)).*log(abs(X)))) .* ((1/N).*sum(abs(X).^(2*a(i-1)))) .* ((2/N).*sum(abs(X).^(a(i-1))))) / (((1/N).*sum(abs(X).^(a(i-1)))).^4) ) - 1;

    % Netwon-Raphson Search.
    a(i) = a(i-1) - Zn(i-1)./ZnP(i-1);
end

estimate = a(end);
total = a;
return;
end
