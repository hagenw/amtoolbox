function outsig = pinknoise(nsamples)
% PINKNOISE Generates a pink noise signal
%   Usage: outsig = pinknoise
%
%   Input parameters:
%       nsamples    - length of the noise (samples)
%
%   Output parameters:
%       outsig      - nsamples x 1 signal vector
%
%   PINKNOISE(samplesn) generates an one channel output signal containing pink 
%   noise (1/f spectrum) with the length of nsamples.
%
%   See also:
%
%R FIXME wikipedia (script after little2007)
%

%   AUTHOR: Hagen Wierstorf


% ------ Checking of input parameter -------------------------------------

error(nargchk(1,1,nargin));

if ~isnumeric(nsamples) || ~isscalar(nsamples) || nsamples<=0
    error('%s: nsamples has to be a positive scalar.',upper(mfilename));
end


% ------ Computation -----------------------------------------------------
fmax = floor(n/2)-1;
f = (2:(n2+1))';
% 1/f amplitude factor
a = 1/sqrt(f);
% Random phase
p = randn(fmax,1) + 1i * randn(fmax,1);
sig = a .* p;
% Create the whole signal (needed for ifft)
d = [1; sig; 1/(fmax+2); flipud(conj(sig))];
% IFFT to get the time signal
outsig = real(ifft(d));
% Scale output
outsig = outsig ./ (max(abs(outsig(:)))+eps);
