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
[fine_ipd,mod_ipd,fine_itd,mod_itd,ild,fine_ic,mod_ic,fine_f_inst,mod_f_inst,cfreqs] = dietz(signal,fs);

% convert interaural information into azimuth
angl=itd2angle(fine_itd,ild(:,1:12),fine_f_inst,lookup,2.5);

% plot
figure;hist(angl(fine_ic>ic_threshold&[diff(fine_ic)>0; zeros(1,12)]),91);
t=(1:length(signal))*1/fs;
figure;plot(t,angl(:,cn));
hold on;
plot(t(fine_ic(:,cn)>ic_threshold),angl(fine_ic(:,cn)>ic_threshold,cn),'r.');

