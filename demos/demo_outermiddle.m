%DEMO_OUTERMIDDLE  Outer and middle ear filters.
%
%  This demo shows the frequency response of the outer and middle ear
%  filters.
%
%  FIGURE 1 Outer ear
%
%    This figure shows the frequency response of a filter modelling the
%    combined effect of the headphones and the outer ear canal.
%
%  FIGURE 2 Middle ear
%
%    This figure shows the frequency response of a filter modelling the
%    effect of the middle ear.
%
%  See also: headphonefilter, middleearfilter

% Sampling frequency to use. Must be higher than 20000 to show the full
% range of the data defining the filters.
fs=22050;

% Calculate the filters.
bout = headphonefilter(fs);
bmid = middleearfilter(fs);

% Manually calculate the frequency response
fout = 20*log10(abs(fftreal(bout)));

% Half the filter length.
n2=length(fout);

% x-values for plotting.
xplot=linspace(0,fs/2,n2);

figure(1);
semiaudplot(xplot,fout);
xlabel('Frequency in Hz (on erb scale)');
ylabel('Magnitude response (dB)');
title('Magnitude response of headphone+outer ear filter.');

figure(2);
% Manually calculate the frequency response
fmid = 20*log10(abs(fftreal(bmid)));

semiaudplot(xplot,fmid);
xlabel('Frequency in Hz (on erb scale)');
ylabel('Magnitude response (dB)');
title('Magnitude response of middle ear filter.');