% dietz_demo.m
%
% loads a simulated binaural recording of multiple persons talking and
% calculates and plots azimuth histograms and the azimuth over time
% for the 1-kHz channel (as in Dietz 2011 paper)

load speakers24; %two talkers
%load speakers123; % 3 talkers
%load speakers1234; % 4 talkers
%load speakers12345; % five talkers

fs=16000;
load azi_lookup_4;
ic_threshold=0.98;
cn = 10; % channel number for time plot

% run dietz model on signal
[hairc_fine, hairc_mod, fc, hairc_ild]=dietz2011(signal,fs);

% convert interaural information into azimuth
angl=itd2angle(hairc_fine.itd_lp,hairc_ild(:,1:12),hairc_fine.f_inst,lookup,2.5);

% plot
figure;hist(angl(hairc_fine.ic>ic_threshold&[diff(hairc_fine.ic)>0; zeros(1,12)]),91);
t=(1:length(signal))*1/fs;
figure;plot(t,angl(:,cn));
hold on;
plot(t(hairc_fine.ic(:,cn)>ic_threshold),angl(hairc_fine.ic(:,cn)>ic_threshold,cn),'r.');
