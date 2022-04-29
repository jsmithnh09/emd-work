function y = emdK(IMF)
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
% See Also: EMDB, GGDKURT
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
% Perform the actual IMF filtering.
%**************************************************************************

% if a single vector, apply EMD to extract IMF's.
if (isvector(IMF))
    IMF = emd_rilling(IMF);
end
[m, ~] = size(IMF);

valid = [false; true(m-1, 1)];
for imfInd = 2:m
    % determine 4th moment estimate.
    k = ggdcsk(IMF(imfInd, :));
    
    % linear lookup.
    [~, tableInd] = min(abs(pKurt.K - k));
    b_est = pKurt.B(tableInd);
    valid(imfInd) = (b_est < 1) || (b_est > 3);
end

% logical indexing to sum only the "valid" IMFs outside (1, 3) region.
y = sum(IMF(valid, :), 1);

end % emdK

