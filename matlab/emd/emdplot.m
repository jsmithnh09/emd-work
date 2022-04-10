function varargout = emdplot(imfmat, inds)
%EMDPLOT plots the IMFs present in the IMF matrix.
%
%   EMDPLOT(IMFMAT) will generate a plot for each of the IMFs in a single
%       figure.
%
%   EMDPLOT(IMFMAT, NUMS) will only plot the IMFs specified in the array
%       NUMS to simplify observations.
%
%   HFIG = EMDPLOT(...) returns the handle to the primary GUI figure, in
%       case the contents of the figure need to be modified.
%   
% See Also: EMD_RILLING

[M, N] = size(imfmat);
tmprange = 1:M;
if (nargin > 1)
    if (any(~ismember(inds, tmprange)))
        error('Specified IMF indices are not valid.');
    end
    indrange = inds;
else
    indrange = tmprange;
end
nIMFs = length(indrange);

H = figure('Visible', 'off');
t = 1:N;
for idx = 1:length(indrange)
    subplot(nIMFs, 1, idx);
    plot(t, imfmat(indrange(idx), :)); grid('on');
    xlim([0, t(end)]); 
    ylabel(sprintf('IMF %d', indrange(idx)));
end
set(H, 'Visible', 'on');
if (nargout >= 1)
    varargout{1} = H;
end

end % emdplot