function outsig = transposedtone(siglen,fc,fm,fs,varargin)
%TRANSPOSEDTONE  Transposed tone test stimuli
%   Usage:  ts = transposedtone(fc,fm,dur,fs);
%
%   Input parameters:
%     siglen  : Length of signal
%     fc      : Vector of carrier frequencies (Hz)
%     fm      : Vector of modulation frequencies (Hz)
%     fs      : Sampling frequency (Hz)
% 
%   Output parameters:
%     outsig  : transposed tone (column vector)
%
%   TRANSPOSEDTONE(siglen,fc,fm,dur,fs) generates a transposed tone test
%   stimuli as defined in Kolrausch et. al (1997).
%
%   By default, the output is normalized to have an RMS value of 1, but this
%   can be changed by passing any of the flags from the NORMALIZE function.
% 
%   Some example parameters as used in a study by Santurette
%   
%C     siglen = 44100;
%C     fc     = 5000;
%C     fm     = 435;
%C     fs     = 44100;
%C     outsig = transposedtone(fc,fm,dur,fs);
%
%   References:kohlrausch1997detection oxenham2004correct

% Author: Sébastien Santurette  2009

if nargin<4
  error('Too few input parameters.');
end;

definput.import={'normalize'};
definput.importdefaults={'rms'};
[flags,keyvals]=ltfatarghelper({},definput,varargin);


% time vector
t = (0:siglen-1)/fs;

% number of tones
n = length(fc);
                       
outsig = zeros(siglen,1);
for ii = 1:n
  % carrier tone
  carrier = cos(2*pi*fc(ii)*t);        

  % "modulator"   
  hrsine = sin(2*pi*fm(ii)*t);         

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

outsig=normalize(outsig,flags.norm);

%OLDFORMAT
