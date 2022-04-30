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
% Contact: Jordan R. Smith

if (~ischar(type))
  error('common:acn:InvalidInput', ...
    'Input TYPE must be type CHAR.');
elseif (~isscalar(SNRdB) || ~isnumeric(SNRdB))
  error('common:acn:InvalidInput', ...
    'Input SNRdB must be a scalar power decibel.');
end
SNRlin = 10^(SNRdB/10); % SNR measured as Ps/Pn, or Ps/var(n) if zero-mean.

% re-orient input signal.
x = x(:);
if (~isvector(x))
    error('common:acn:InvalidInput', ...
        'Input SIGNAL must be a unidimensional signal.');
end
npts = numel(x);

% generate the noise.
switch type
    case 'white'
        n = randn(npts,1);
    case 'pink'
        n = pinknoise(npts);
    case {'violet', 'purple'}
        n = violetnoise(npts);
    case 'red'
        n = rednoise(npts);
    case 'blue'
        n = bluenoise(npts);
    otherwise
        error('common:acn:InvalidInput', ...
            'Unrecognized noise type "%s".', type);
end % switch

% unity variance, zero-mean.
n = n(:);
n = (n - mean(n))/std(n, 1);

% delegate the appropriate noise level. Ps/Pn drops (1/N) term.
Ps = sum(x.^2);
Pn = sum(n.^2);
alpha = sqrt(Ps / (SNRlin * Pn));
y = x + alpha * n;
return;
end % acn

%**************************************************************************
%               rednoise
%**************************************************************************
function y = rednoise(N)

% function: y = rednoise(N)
% N - number of samples to be returned in a row vector.
% Y - a row vector of red (Brownian) noise samples.
%
% the function generates a sequence of red (Brownian) noise samples.
% In terms of power at a constant bandwidth, red noise falls off at 6 dB/oct, i.e. 20 dB/dec.

% define the length of the vector and ensure that the M is even.
if rem(N, 2)
    M = N+1;
else
    M = N;
end

% generate white noise
x = randn(1, M);

% FFT
X = fft(x);

% prep a vec with freq indexes
nPts = M/2 + 1;
n = 1:nPts;

% manipulate left half of the spectrum so the PSD is propertional to the
% frequency by a factor of 1/(f^2), or the amplitudes are proportional to 1/f
X = X(1:nPts);
X = X./n;

% prep the right half of the spectrum, a conj. copy of the left
% except the DC component and the Nyquist component, (0 < x < Nq)
X = [X conj(X(end-1:-1:2))];

% IFFT
y = real(ifft(X));

% ensure that the length of y is N.
y = y(1, 1:N);

end % rednoise

%**************************************************************************
%                           bluenoise
%**************************************************************************
function y = bluenoise(N)

% function: y = bluenoise(N) 
% N - number of samples to be returned in a row vector
% y - a row vector of blue noise samples

% The function generates a sequence of blue noise samples. 
% In terms of power at a constant bandwidth, blue noise increase in at 3 dB/oct, i.e. 10 dB/dec. 

% difine the length of the vector and
% ensure that the M is even, this will simplify the processing
if rem(N, 2)
    M = N+1;
else
    M = N;
end

% generate white noise 
x = randn(1, M);

% FFT
X = fft(x);

% prepare a vector with frequency indexes 
nPts = M/2 + 1;             % number of the unique fft points
n = 1:nPts;                 % vector with frequency indexes 

% manipulate the left half of the spectrum so the PSD
% is proportional to the frequency by a factor of 1/f, 
% i.e. the amplitudes are proportional to 1/sqrt(f)
X = X(1:nPts);  
X = X.*sqrt(n);

% prepare the right half of the spectrum - a conjugate copy of the left one,
% except the DC component and the Nyquist component - they are unique
% and reconstruct the whole spectrum
X = [X conj(X(end-1:-1:2))];

% IFFT
y = real(ifft(X));

% ensure that the length of y is N
y = y(1, 1:N);

end % bluenoise

%**************************************************************************
%               violetnoise
%**************************************************************************
function y = violetnoise(N)

% function: y = violetnoise(N) 
% N - number of samples to be returned in a row vector
% y - a row vector of violet noise samples

% The function generates a sequence of violet (purple) noise samples. 
% In terms of power at a constant bandwidth, violet noise increase in at 6 dB/oct, i.e. 20 dB/dec. 

% difine the length of the vector and
% ensure that the M is even, this will simplify the processing
if rem(N, 2)
    M = N+1;
else
    M = N;
end

% generate white noise 
x = randn(1, M);

% FFT
X = fft(x);

% prepare a vector with frequency indexes 
nPts = M/2 + 1;     % number of the unique fft points
n = 1:nPts;         % vector with frequency indexes 

% manipulate the left half of the spectrum so the PSD
% is proportional to the frequency by a factor of f^2, 
% i.e. the amplitudes are proportional to f
X = X(1:nPts);
X = X.*n;

% prepare the right half of the spectrum - a conjugate copy of the left one,
% except the DC component and the Nyquist component - they are unique
% and reconstruct the whole spectrum
X = [X conj(X(end-1:-1:2))];

% IFFT
y = real(ifft(X));

% ensure that the length of y is N
y = y(1, 1:N);


end % violetnoise
%**************************************************************************
%               pinknoise
%**************************************************************************
function y = pinknoise(N)

% function: y = pinknoise(N) 
% N - number of samples to be returned in a row vector
% y - a row vector of pink (flicker) noise samples

% The function generates a sequence of pink (flicker) noise samples. 
% In terms of power at a constant bandwidth, pink noise falls off at 3 dB/oct, i.e. 10 dB/dec. 

% define the length of the vector and
% ensure that the M is even, this will simplify the processing
if rem(N, 2)
    M = N+1;
else
    M = N;
end

% generate white noise
x = randn(1, M);

% FFT
X = fft(x);

% prepare a vector with frequency indexes 
nPts = M/2 + 1;     % number of the unique fft points
n = 1:nPts;         % vector with frequency indexes 

% manipulate the left half of the spectrum so the PSD
% is proportional to the frequency by a factor of 1/f, 
% i.e. the amplitudes are proportional to 1/sqrt(f)
X = X(1:nPts);      
X = X./sqrt(n);

% prepare the right half of the spectrum - a conjugate copy of the left one,
% except the DC component and the Nyquist component - they are unique
% and reconstruct the whole spectrum
X = [X conj(X(end-1:-1:2))];

% IFFT
y = real(ifft(X));

% ensure that the length of y is N
y = y(1, 1:N);

end % pinknoise
