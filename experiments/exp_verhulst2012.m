function exp_verhulst2012
%EXP_VERHULST2012 Compute figurs from the Verhulst paper
%
%    this script compute and plot figures 2.a and 2.c of: Verhulst, Sarah,
%   Torsten Dau, and Christopher A. Shera.  "Nonlinear time-domain cochlear
%   model for transient stimulation and human otoacoustic emission."  The
%   Journal of the Acoustical Society of America 132.6 (2012): 3842-3848.
%

%   AUTHOR: Alessandro Altoè 

% close all;
% clear all;
% clc;
fs=48000;
f0=1e3;
t=0:1/fs:0.05;
levels_n0=9;

%% plot fig. 2c
%ochlear excitation patterns calculated as the rms level of displacement per cochlear section,
%for stimulation with a pure tone of 1 kHz with stimulus intensities between 10 and 90 dB SPL
sig=ones(levels_n0,1)*(sin(2*pi*f0*t));
onset_duration=1e-2*fs;
sig(:,1:round(onset_duration)+1)=sig(:,1:round(onset_duration)+1).*(ones(levels_n0,1)*(0:1.0/round(onset_duration):1)); %fade in 10ms
spl=10:80/(levels_n0-1):90;
fc='all';
norm_Rms=ones(levels_n0,1);
irr=ones(levels_n0,1)';
[v,y,e,cf]=verhulst2012(sig,fs,fc,spl,norm_Rms,1,irr);
Y_rms=mag2db(rms(y(:,:,:)./0.01));
leg_str=repmat('dB',levels_n0,1);
figure;
for i=1:levels_n0
semilogx(cf(2:length(cf)),(Y_rms(:,2:length(cf),i)));
hold all;
end
legend((horzcat(int2str(spl'),leg_str)));
grid on;
set(gca, 'xdir','reverse');
set(gca,'Xtick',[250 500 1000 2000 4000 8000 16000]);
set(gca,'XTickLabel',{'0,25','0.5','1','2','4','8','16'});
xlabel('Center frequency (KHz)');
ylabel('Y_{rms} (dB re .01)');
axis([250 8000 -250 -100]);
title('BM displacement in response to 1 KHz sinsusoid');

%% plot fig 2a
%BM displacement simulated for the 1-kHz cochlear CF location for clicks with intensities between 0 and 90 dB peSPL. The displacements were normalized 
%by the pressure at the stapes of the cochlea such that compression is observed as a reduction of the IR amplitude 
clear fc v y e cf fc
fc=[1e3];
norm_Rms=norm_Rms.*0;
sig=sig.*0;
sig(:,8:9)=2; %impulse (the value 2 is to have 2 peak-to-peak value)
[v,y,e,cf]=verhulst2012(sig,fs,fc,spl,norm_Rms);
figure;
for i=1:levels_n0
plot(t*1e3,y(:,1,i)./max(abs(y(:,i))));
hold all;
end
legend((horzcat(int2str(spl'),leg_str)));
grid off;
xlabel('Time (ms)');
ylabel('Normalized y_{BM}');
axis 'tight';
title('Displacement in response to a click');
xlim([0 25]);
