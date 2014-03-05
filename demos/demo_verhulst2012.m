%DEMO_VERHULST2012 Demo of the OtoacousticEmission
%  
%   This script computes and plot the otoacoustic emission for a 500 Hz sinusoid,
%   This is generated using the cochlear model described in verhulst2012
%   In particular otoacoustic emissions are computed as the signal difference between the
%   sound pressure at the middle ear with model non linearities and irregularities enabled,
%   and the sound pressure at the middle ear in case of linear model.
%
%   XXX Descriptions of the figures are missing

fs=48000;
t=0:1.0/fs:0.05;
insig=zeros(length(t),2);
f0=500;
insig(:,1)=sin(2*pi*f0*t);
insig(:,2)=insig(:,1);
spl=[60,60];
irron=[1,0];
normalizeRms=[1 1];
subjectNo=rand();
fc='all';

[V,Y,E,CF]=verhulst2012(insig,fs,fc,spl,normalizeRms,subjectNo,irron);

OtoacousticEmissionPa=E(:,1)-E(:,2);

figure; 
plot(t.*1e3,OtoacousticEmissionPa);
grid on;
xlabel('time (ms)');
ylabel('Emission (Pa)');
title('Otoacoustic Emission for a 500 Hz sinsuoid');