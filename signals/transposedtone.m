function outsig = transposedtone(fc,fm,dur,ear,tslevel)
%TRANSPOSEDTONE  Transposed tone test stimuli
%   Usage:  ts = transposedtone(fc,fm,dur,ear,tslevel);
%
%   Input parameters:
%     fc      : Vector of carrier frequencies (Hz)
%     fs      : Sampling frequency
%     fm      : Vector of modulation frequencies (Hz)
%     dur     : Stimulus duration (s)
% 
%   Output parameters:
%     outsig  : transposed tone (column vector)
%
%R  kohlrausch1997detection oxenham2004correct

% Author: Sébastien Santurette  2009

% time vector
t = (0:dur-1)/fs;

% number of samples
ns = dur*fs;
                       
% number of tones
n = length(fc);
                       
outsig = zeros(ns,1);
for ii = 1:n
  % carrier tone
  carrier = cos(2*pi*fc(i)*t);        

  % "modulator"   
  hrsine = sin(2*pi*fm(i)*t);         

  % half-wave rectify the modulator
  hrsine = max(0,hrsine);             
  
  % Compute coefficients for low-pass filtering
  fcuth = .2*fc(ii)/(fs/2);
  [b,a] = butter(4,fcuth,'low');

  % Low-pass filter the modulator
  hrsine = filter(b,a,hrsine);

  % Add the carrier*modulator to the transposed tone
  outsig = outsig+(carrier.*hrsine)';
end;

% Normalize the output
outsig = outsig/rms(outsig);