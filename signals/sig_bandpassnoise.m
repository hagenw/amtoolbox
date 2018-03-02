function outsig = sig_bandpassnoise(fc,fs,dur,l,bw)
%sig_bandpassnoise  Generates a bandpass-noise-type masker
%
%   Usage:   outsig = sig_bandpassnoise(fc,fs,dur,l,bw)
%
%   `outsig = sig_bandpassnoise(fc,fs,dur,l,bw)` generates a bandpass-noise
%   masker with duration `dur` (in sec), overall level `l` (in dB SPL) and
%   with a sampling rate of `fs` Hz at the center frequency fc. The 
%   bandpass has a bandwidth of `bw` in Hz. The generation is in the
%   frequency domain.
%
%   Input parameters:
%       fc:     center frequency of the bandpass
%       fs:     sampling frequency
%       dur:    duration of the bandpassnoise in s
%       l:      overall level in (db SPL) with a dB-offset of 0 dB
%       bw:     bandwidth of the noiseband
%
%   Output parameters:
%       outsig: the generated bandpassnoise
%
%   Author:     Martina Kreuzbichler

if nargin<5
  error('%s: Too few input parameters.',upper(mfilename));
end;


% Make sure the length is even
n = round(dur*fs/2)*2;

% Caluculation of variables
fl = fc-bw/2;
fh = fc+bw/2;
resolution = fs/n;
flpin = round(fl/resolution);
fhpin = round(fh/resolution);

%special case when bandwidth is twice the centerfrequency
if flpin == 0
    flpin = 1;
end




%% Generate broadband Gaussian noise

% Make sure the length is even
n = round(dur*fs/2)*2;

%% Generation in the frequency domain
a = zeros(n/2+1,1);
a_dB = a;
inds = flpin:fhpin;
a_dB(inds) = 1;
a(inds) = 10.^(a_dB(inds)/20)*sqrt(n)/2;
p = randn(n/2+1,1) + 1i*randn(n/2+1,1);
sig = bsxfun(@times,a,p);
noisebp = ifftreal(sig,n);
outsig = setdbspl(noisebp,l,'dboffset',0);  
 

%% Plot results (uncomment if needed)

% fax = 0:resolution:fs-resolution;
% plot(fax, abs(fft(noisebplv)))

