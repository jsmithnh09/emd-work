function y = acn(x, SNRdB, type)
% Additive Colored Noise with user specified dB contamination.
%
% Y = ACN(SIGNAL, SNR, TYPE)
%   SIGNAL (double) is the input signal to contaminate.
%   SNR (double) is a scalar of power dB contamination.
%   TYPE (char) is the type of noise to apply. Valid inputs include:
%       'WHITE'     -- standard normally distributed white noise.
%       'RED'       -- brownian noise, (random motion of particles.)
%       'BLUE'      -- blue noise found in Cherenkov Radiation
%       'VIOLET'    -- or purple, same as acoustic thermal noise in water.
%       'PINK'      -- reference noise in audio engineering.
%
% Spectral manipulation based off of alpha scaling described with [1].
%
% References:
% [1] H. Zhivomirov. A Method for Colored Noise Generation. Romanian Journal 
%     of Acoustics and Vibration, ISSN: 1584-7284, Vol. XV, No. 1, pp. 14-19, 2018. 
%     (http://rjav.sra.ro/index.php/rjav/article/view/40/29)
%
% Contact: Jordan R. Smith

if (~ischar(type))
  error('common:acn:InvalidInput', ...
    'Input TYPE must be type CHAR.');
elseif (~isscalar(SNRdB) || ~isnumeric(SNRdB))
  error('common:acn:InvalidInput', ...
    'Input SNRdB must be a scalar power decibel.');
end
SNRlin = db2pow(SNRdB);

% re-orient input signal.
x = x(:);
if (~isvector(x))
    error('common:acn:InvalidInput', ...
        'Input SIGNAL must be a unidimensional signal.');
end
N = numel(x);

% generate the noise.
switch type
    case 'white'
        alpha = 0;
    case 'pink'
        alpha = -1;
    case {'violet', 'purple'}
        alpha = 2;
    case 'red'
        alpha = -2;
    case 'blue'
        alpha = 1;
    otherwise
        error('common:acn:InvalidInput', ...
            'Unrecognized noise type "%s".', type);
end % switch

% convert PSD to Amplitude Density (div.by 2.)
alpha = alpha / 2;
w = randn(N,1);

% manipulate spectral roll-off in the left-sided spectrum. no need to scale
% since we'll mirror before going back to time.
if (alpha ~= 0)
    M = ceil((N/2)+1);      % FFT point length
    W = fft(w);             % get the frequency content
    fvec = (1:M).';         % freq. indexes
    W = W(fvec);            % single-sided spectrum
    W = W.*(fvec.^alpha);   % spectral amp. proportional to factor f^alpha

    % even-length FFT, include Nyquist point.
    if (rem(N, 2))
        W = [W; conj(W(end:-1:2))];
    else
        % odd length FFT, exclude Nyquist.
        W = [W; conj(W(end-1:-1:2))];
    end
    n = real(ifft(W));
else
    % white-noise case.
    n = w;
end

% zero mean, unity variance.
n = n - mean(n);
n = n / std(n);
n = n(:);

% delegate the appropriate noise level. Ps/Pn drops (1/N) term.
Ps = sum(x.^2);
Pn = sum(n.^2);
beta = sqrt(Ps / (SNRlin * Pn));
y = x + (beta * n);
return;

end % acn
