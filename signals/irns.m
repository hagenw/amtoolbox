function outsig = irns(delay,iterations,gn,siglen,fs);
%IRNS	Generate iterated rippled noise
% 	Usage: outsig=irns(delay,iterations,gn,siglen,fs)
%
%   IRNS(delay,iterations,gn,siglen,fs) generates a signal consisting of white 
%   noise, which is added iterations times to himself with a delay of delay (in 
%   ms). 
%
%   Options:
%  		delay      	delay in ms of the time-shifted noise adding process
%  		iterations 	number of iterations of the adding process
%  		gn       		relative gain of irn
%  		siglen	   	signal length of irn in samples
%  		fs         	sampling rate in Hz
%
% Try
% fs = 44100;
% signal = irns(4,6,1,fs,fs);
% sound(signal,fs)
%
% See also: irno

% AUTHOR: Hagen Wierstorf, Daniel Pressnitzer, Stefan Uppenkamp

% ------ Checking of input parameters ---------

error(nargchk(5,5,nargin));

% ------ Computation --------------------------

% Frequency to which the delay corresponds
freq = 1000/delay;
% Number of samples for the noise (slightly longer to avoid circular iterations 
% (S. Uppenkamp)
noiselen = siglen + round(iterations*fs/freq);
% Number of samples for the delay
delaylen = round(fs/freq);
% Create white noise
noisesig = randn(1,noiselen);

% Iterate delay and add n times
for ii = 1:iterations
  dnoise(delaylen+1:noiselen) = noisesig(1:noiselen-delaylen);
  dnoise(1:delaylen) = noisesig(noiselen-delaylen+1:noiselen);
  noisesig = noisesig + gn*dnoise;
end

% Take first bit of result as IRN
outsig = noisesig(1:siglen);
% Scale to RMS of 1
outsig = outsig/rms(outsig);

