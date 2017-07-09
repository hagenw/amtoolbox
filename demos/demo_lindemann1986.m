%DEMO_LINDEMANN1986 Demo of the Lindemann binaural model
%
%   This script generate a figure showing the result of the lindemann
%   binaural model for a 2 Hz binaural modulated sinusoid with a frequency of
%   500 Hz.
%
%   .. figure::
%
%      Binaural modulated sinusoid
%
%      This figure shows the binaural activity map for one frequency channel of
%      the lindemann binaural model for a sinusoid with a binaural modulation
%      rate of 2 Hz.
%
%   .. figure::
%
%      Sinusoid with ITD
%
%      This figure shows the result of the Lindemann 1986 binaural model
%      averaged over time for the desired frequency channel for a sinusoid
%      with an ITD of 0.3 ms.
%
%   See also: lindemann1986, lindemann1986_bincorr, plot_lindemann1986


% Sampling rate
fs = 44100;
% Frequency of the sinusoid
f = 500;


% ------ Fig 1. ----------------------------------------------------------

figure(1);

% Binaural modulation frequency
mf = 2;
% Generate 1~s binaural modulated sinusoid
sig = sig_lindemann1986(f,mf,fs);

% Model paramter (Note: T_int (ms) should be a multiple of 1000/f == 2)
% Begin of the storage of the cross-correlation is set to 1, because we have a
% non-stationary signal

% Calculate binaural cross-correlation
[cc,t] = lindemann1986(sig,fs,'T_int',6);

% Set title string for the plot
tstr = sprintf(['Binaural modulated sinusoid\nf = %i Hz\nf_m = %i Hz\n',...
    'fc = %i\n'],f,mf,round(freqtoerb(f)));
% Plot frequency channel 11, due to round(freqtoerb(500))==11
plot_lindemann1986(cc,t,'fc',f,'title',tstr);


% ------ Fig 2. ----------------------------------------------------------

figure(2);

% Generate an sinusoid with a ITD
itd = 0.3; % (ms)
sig = sig_itdsin(f,itd,fs);
sig = sig(1:fs/2,:);

% Calculate binaural cross-correlation using the 'stationary' mode and
% Lindemanns default model parameter
[cc,t] = lindemann1986(sig,fs,'stationary');

% Set title string for the plot
tstr = sprintf('Sinusoid with an ITD\nf = %i Hz\nitd = %.1f ms\nfc = %i\n',...
    f,itd,round(freqtoerb(f)));
% Plot frequency channel 11, due to round(freqtoerb(500))==11
plot_lindemann1986(cc,t,'fc',f,'title',tstr);
