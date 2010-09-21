% This program illustrates/test the function modspecgram

% Parameters used below
fm = 50;    % Modulation frequency in (1)
l = 2;      % Length of the signal in (1)in seconds
fs = 44100; % Sampling frequency in (1)
fc = 5000;  % Carrier frequency in (3)
M = 500;
a = 30;
g = 'gauss';

%% 1) Modulation spectrogram of a modulated noise

% Time vector
t = 0:1/fs:l;
n = length(t);

% Noise and modulated noise
noise = 1-2*randn(1,n);
%noise =cos(2*pi*2000*t);
modnoise = noise.*(1+cos(2*pi*t*fm));

figure(1)
modspecgram(modnoise,fs,90)
title('Sinusoidaly Modulated Noise')


%% 2) Modulation spectrogram of two speech samples

figure(2)
modspecgram(greasy,16000,60,500)
title('Greasy')

figure(3)
modspecgram(greasylong,8000,60,100)
title('Greasy (long)')

%% 3) Sinusoidal carrier

s = sin(2*pi*t*fc);
smod = s.*(1+0.5*cos(2*pi*t*fm));

figure(4)
modspecgram(s,fs,50,2*fm,'fmax',2*fc)
