% ex_verhulst
% this script compute and plot figures 2.a and 2.c of verhulst2012
%
close all;
clear all;
clc;
fs=96000;
f0=1000;
t=0:1/fs:0.02;
levels_n0=9;

%% plot fig. 2c
%ochlear excitation patterns calculated as the rms level of displacement per cochlear section,
%for stimulation with a pure tone of 1 kHz with stimulus intensities between 10 and 90 dB SPL
sig=ones(levels_n0,1)*(sin(2*pi*f0*t));
sig(:,1:round(fs*4e-3)+1)=sig(:,1:round(fs*4e-3)+1).*(ones(levels_n0,1)*(0:1.0/round(fs*4e-3):1)); %fade in
spl=10:80/(levels_n0-1):90;
fc='all';
norm_Rms=ones(levels_n0,1);
irr=ones(levels_n0,1)';
[v,y,e,cf]=verhulst2012(sig,fs,fc,spl,norm_Rms,1,irr);
Y_rms=mag2db(rms(y(:,:,:)./0.1));
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
ylabel('Y_{rms} (dB re .1)');
axis([250 8000 -250 -100]);
title('BM displacement in response to 1 KHz sinsusoid');

%% plot fig 2a
%BM displacement IRs simulated for the 1-kHz cochlear CF location for clicks with intensities between 0 and 90 dB peSPL. The IRs were normalized 
%by the pressure at the stapes of the cochlea such that compression is observed as a reduction of the IR amplitude 
clear fc v y e cf fc
fc=[1e3];
norm_Rms=norm_Rms.*0;
sig=sig.*0;
sig(:,8:12)=2; %impulse
[v,y,e,cf]=verhulst2012(sig,fs,fc,spl,norm_Rms);
norm_factor=max(max(max(abs(y))));
figure;
for i=1:levels_n0
plot(t*1e3,y(:,1,i)./max(abs(e(:,i))));
hold all;
end
legend((horzcat(int2str(spl'),leg_str)));
grid on;
xlabel('time (ms)');
ylabel('Normalized displacement');
axis 'tight';
title('IRs simulated for the 1-kHz cochlear CF location for clicks with intensities between 0 and 90 dB peSPL');