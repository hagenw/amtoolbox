%DEMO_LINDEMANN
%
%   This script generate a figure showing the result of the lindemann
%   binaural model for a 2 Hz binaural modulated sinusoid.

% Sampling rate
fs = 44100;
% Frequency of the sinusoid
f = 500;
% Binaural modulation frequency
mf = 2;

% Generate binaural modulated sinusoid
sig = bmsin(f,mf,fs);

% Model paramter
w_f = 0.035; c_s = 1;

crosscorr = lindemann(sig,0.035,1,fs);
plotlindemann(crosscorr);