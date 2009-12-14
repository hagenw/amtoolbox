function outsig = bandpassnoise(siglen,flow,fhigh,fs);
%BANDPASSNOISE  Generate band pass noise
%   Usage: outsig=bandpassnoise(siglen,flow,fhigh,fs);
%
%   BANDPASSNOISE(siglen,flow,fhigh,fs) generates a signal consisting of
%   white noise in a frequency band from flow to fhigh. The signal is
%   samples at a sample rate of fs Hz.
%
%   The signal is normalized to an RMS value of 1.

randsig = fft(randn(siglen,1));

% Round them to integer values.
flow = round(flow*siglen/fs);
fhigh = round(fhigh*siglen/fs);

% HACK: if lowpass ( flow = 0) index would be greater than len (len +1)
flow=max(flow,1);

outsig = zeros(siglen,1);
outsig(flow+1:fhigh+1) = randsig(flow+1:fhigh+1);

outsig(siglen-fhigh+1:siglen-flow+1) = randsig(siglen-fhigh+1:siglen-flow+1);

outsig = real(ifft(outsig));

outsig = outsig/rms(outsig);




