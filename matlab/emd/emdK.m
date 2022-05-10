function y = emdK(IMF,varargin)
% EMD-Kurtosis utilizes 4th moment observations to estimate the shape
% parameter.
%
%   Y = EMDK(IMF)
%     IMF (double) is the M-by-N matrix from EMD.
%     Y (double) is an N-by-1 denoised vector.
%
%   [~] = EMDK -CLEAR removes the static memory of the theoretical kurtosis
%       vector. Keeping this in memory saves time of computation when
%       performing Monte-Carlo iterations.
%
%   Y = EMDK(IMF, ... 'LowerBound', 1, 'UpperBound', 3)
%       Accepts additional keyword argumets to modify the GGD distribution
%       boundaries via "LOWERBOUND" and "UPPERBOUND".
%
%   Y = EMDK(IMF, ..., 'BoundsCheck', 0, 'ZScore', 1.96)
%       Accepts additional keyword argument "BOUNDSCHECK" to determine if
%       we want to apply the normality test when we're within +/- tolerance
%       of the GCM confidence interval utilized in EMD-Beta. The "ZSCORE"
%       can be modified as well, using 95% by default.
%
% See Also: EMDB, GGDKURT, EMD_RILLING, OPTPARSE
%
% Contact: Mahdi Al-Badrawi, Jordan Smith

%**************************************************************************
% Store the kurtosis table for lookup.
%**************************************************************************
persistent pKurt;
if (isempty(pKurt))
    Baxis = 0.1:0.0001:4;
    Kcoeff = ggdkurt(Baxis);
    pKurt = struct('B', Baxis, 'K', Kcoeff);
end
if (ischar(IMF) && strcmp(IMF, '-clear'))
    clear('pKurt'); 
    y = [];
    return;
end

%**************************************************************************
% Check keyword arguments.
%**************************************************************************
opt = struct(...
    'BoundsCheck', 0, ...
    'LowerBound', 1, ... 
    'UpperBound', 3, ... 
    'ZScore', 1.96);
opt = optparse(opt, varargin{:});

% if a single vector, apply EMD to extract IMF's.
if (isvector(IMF))
    IMF = emd_rilling(IMF);
end
[m, n] = size(IMF);

% perform the IMF validation.
valid = false(m, 1);
for imfInd = 2:m
    
    % determine 4th moment estimate.
    r = IMF(imfInd, :);
    k = ggdcsk(r);
    
    % linear lookup.
    [~, tableInd] = min(abs(pKurt.K - k));
    b_est = pKurt.B(tableInd);

    % bounds check.
    valid(imfInd) = local_bounds_check(r, b_est, n, opt);

end

% logical indexing to sum only the "valid" IMFs outside (1, 3) region.
y = sum(IMF(valid, :), 1);

end % emdK

%**************************************************************************
% BOUNDS_CHECK
%   Determine if the bounds are reasonable. We have two scenarios we can
%   test:
%       0 - only check if outside the Laplacian and (B > 3) case.
%       1 - get the tolerance from the GCM Confidence Interval and apply
%           Shapiro-Wilk to determine if Normal.
%
% By default, we're using no Normality checks and simply using Kurtosis.
%**************************************************************************
function v = local_bounds_check(samp, shape, N, opts)
    v = false;
    switch opts.BoundsCheck
        case 0
            % no bounds-checking, simply toss result.
            v = (shape < opts.LowerBound) || (shape > opts.UpperBound);
            return;
        case 1
            % confidence interval check based on absolute-value of moments.
            [lbound, ubound] = gcmci(shape, N, opts.ZScore);
            tol = abs(ubound - lbound);
            if ((shape >= opts.LowerBound-tol) && (shape <= opts.UpperBound+tol)) || ...
                    ((shape >= opts.UpperBound-tol) && (shape <= opts.UpperBound+tol))
                [H, ~, ~] = swtest(samp, 0.01, 0);
                v = logical(H);
            end
        otherwise
            error('emdK:bounds_check:UnknownCase:: Unknown shape test-case.');
    end % switch
end