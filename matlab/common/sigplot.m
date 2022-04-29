function varargout = sigplot(varargin)
% SIGPLOT plots signals without the need for wrapper GUI code.
%
%   H = SIGPLOT(X1, X2, ...)
%       X1,2... (vector) are signals of the same length.
%
% See Also:
nsigs = nargin;
for iSig = 1:nargin
    dims = size(varargin{iSig});
    if (numel(dims) ~= 2) || (min(dims) ~= 1)
        error('Each input must be a vector.');
    end
end
h = figure;
for iSig = 1:nsigs
    X = varargin{iSig};
    t = 1:length(X);
    subplot(nsigs, 1, iSig);
    plot(t, varargin{iSig}, 'g'); grid('on');
    set(gca, 'Color', 'k');
    xlim([1 t(end)]);
end
if (nargout == 1)
    varargout{1} = h;
end

end % sigplot