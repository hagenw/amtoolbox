function outsig = viemeister79(insig,fs)
%VIEMEISTER79  The Viemeister (1979) leaky-integrator model
%   Usage: outsig=viemeister79(insig,fs);
%
%   This model is included mostly as a test, as it is so simple.
 
% ------ Checking of input parameters ---------
  
error(nargchk(2,2,nargin));

  
% 4-6 kHz four-pole Butterworth Bandpass (missing)
  
% halfwave rectification
insig = max(insig,0);

% first-order lowpass filter @ 65 Hz
[lp_b,lp_a] = butter(1,65/(fs/2));
insig = filter(lp_b, lp_a, insig);

% ac-coupled rms = std
outsig = std(insig,1);
