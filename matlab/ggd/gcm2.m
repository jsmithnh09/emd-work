function [est, total] = gcm2(x, initial, trials)
% searches GCM function with Newton-Raphson iterative search [1].
%
%   [ESTIMATE, TOTAL] = GCM2(X[, INITIAL, TRIALS])
%       X (vector) is a random sample.
%       INITIAL (scalar) is the first guess of the distribution shape.
%       (Default: m1/sqrt(m2), "m" denotes statistical moments)
%       TRIALS (scalar) are the number of iterations to go through.
%       (Default: 100)
%       TOTAL (scalar) indicates the total # of trials to converge.
%       ESTIMATE (scalar) is the estimated shape parameter.
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
switch nargin
    case 1
        ntrials = 100;
        % offset since we need bias (B >= 2)
        bprev = mean(abs(x)) / std(x) + 3;
    case 2
        ntrials = 100;
        bprev = initial;
    otherwise
        ntrials = trials;
        bprev = initial;
end % switch

% convergence check.
issmall = @(x1,x0)abs(x1-x0) < sqrt(eps);

% compute one step to start loop.
bcur = bprev - (Z(x, bprev) / Zp(x, bprev));
iter = 1;
while(iter < ntrials && ~issmall(bcur, bprev))
    bprev = bcur;
    bcur = bprev - (Z(x, bprev) / Zp(x, bprev));
    iter = iter + 1;
end

% return the estimate and total trials.
total = iter-1;
est = bcur;

end

function y = Z(x, B)
% Z(B) convex equation.
    n = numel(x);
    B2 = 2*B;
    num = (1/n) * sum(abs(x).^B2);
    den = ((1/n) * sum(abs(x).^B))^2;
    y = (num / den) - (B + 1);
end

function y = Zp(x, B)
% Z'(B) prime first-order derivative of convex equation.
    n = numel(x);
    B2 = 2*B;
    ab1 = sum(abs(x).^B);
    ab2 = sum(abs(x).^(B2));

    logb2 = sum( (abs(x).^B2) .* log(abs(x)) );
    logb1 = sum( (abs(x).^B) .*  log(abs(x)) );

    num1 = ((2/n) * logb2) * ((1/n) * ab1)^2;
    num2 = ((1/n) * logb1)*((1/n) * ab2)*((2/n) * ab1);
    den = ((1/n) * ab1)^4;
    

    y = (num1 / den) - (num2 / den) - 1;
end



