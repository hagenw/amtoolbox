%DEMO_LINDEMANN
%
%   This script generate a figure showing the result of the lindemann
%   binaural model for a 2 Hz binaural modulated sinusoid with a frequency of
%   500 Hz.

% Sampling rate
fs = 44100;
% Frequency of the sinusoid
f = 500;


% ------ Fig 1. ----------------------------------------------------------

% Binaural modulation frequency
mf = 2;
% Generate 1~s binaural modulated sinusoid
sig = bmsin(f,mf,fs);

% Model paramter (Note: T_int (ms) should be a multiple of 1000/f == 2)
c_s=0.3; w_f=0.035; M_f=6; T_int=6;
% Begin of the storage of the cross-correlation is set to 1, because we have a
% non-stationary signal
% NOTE: you can also use the 'dynamic' option to achieve this
N_1=1;

% Calculate binaural cross-correlation
[cc,t] = lindemann(sig,fs,c_s,w_f,M_f,T_int,N_1);

% Set title string for the plot
tstr = sprintf(['Binaural modulated sinusoid\nf = %i Hz\nf_m = %i Hz\n',...
    'fc = %i\n'],f,mf,round(freqtoerb(f)));
% Plot frequency channel 11, due to round(freqtoerb(500))==11
plotlindemann(cc,t,tstr,f);


% ------ Fig 2. ----------------------------------------------------------

% Generate an sinusoid with a ITD
itd = 0.3; % (ms)
sig = itdsin(f,itd,fs);
sig = sig(1:fs/2,:);

% Calculate binaural cross-correlation using the 'stationary' mode and
% Lindemanns default model parameter
[cc,t] = lindemann(sig,fs,'stationary');

% Set title string for the plot
tstr = sprintf('Sinusoid with an ITD\nf = %i Hz\nitd = %.1f ms\nfc = %i\n',...
    f,itd,round(freqtoerb(f)));
% Plot frequency channel 11, due to round(freqtoerb(500))==11
plotlindemann(cc,t,tstr,f);
