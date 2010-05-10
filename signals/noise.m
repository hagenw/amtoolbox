function outsig = noise(nsamples)
% PINKNOISE Generates a white noise signal
%   Usage: outsig = noise(nsamples)
%
%   Input parameters:
%       nsamples    - length of the noise (samples)
%
%   Output parameters:
%       outsig      - nsamples x 1 signal vector
%
%   NOISE(samplesn) generates an one channel output signal containing white 
%   noise with the length of nsamples.
%
%   See also:
%
%R FIXME wikipedia
%

%   AUTHOR: Hagen Wierstorf


% ------ Checking of input parameter -------------------------------------

error(nargchk(1,1,nargin));

if ~isnumeric(nsamples) || ~isscalar(nsamples) || nsamples<=0
    error('%s: nsamples has to be a positive scalar.',upper(mfilename));
end


% ------ Computation -----------------------------------------------------
fmax = floor(nsamples/2)-1;
% Amplitude factor
a = 1;
% Random phase
p = randn(fmax,1) + 1i * randn(fmax,1);
sig = a .* p;
% Create the whole signal (needed for ifft)
d = [1; sig; 1/(fmax+2); flipud(conj(sig))];
% IFFT to get the time signal
outsig = real(ifft(d));
% Scale output
outsig = outsig ./ (max(abs(outsig(:)))+eps);
