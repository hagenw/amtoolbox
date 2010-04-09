%DEMO_LINDEMANN
%
%   This script generate a figure showing the result of the lindemann
%   binaural model for a 2 Hz binaural modulated sinusoid with a frequency of
%   500 Hz.

% Sampling rate
fs = 44100;
% Frequency of the sinusoid
f = 500;
% Binaural modulation frequency
mf = 2;

% Generate 1~s binaural modulated sinusoid
sig = bmsin(f,mf,fs);
t = 1;

% Model paramter
c_s=1; w_f=0.035; M_f=6; T_int=5;

% Calculate binaural cross-correlation
crosscorr = lindemann(sig,fs,c_s,w_f,M_f,T_int);

% Set title string for the plot
tstr = sprintf('f = 500 Hz\nf_m = 2 Hz\nfc = 11\n');
% Plot frequency channel 11, due to round(freqtoerb(500))==11
plotlindemann(crosscorr,t,tstr,500);
