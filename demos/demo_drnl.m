%DEMO_DRNL  Widening of filters in the DRNL
%
%   This script displays three auditory spectrograms generated using the
%   DRNL. The input signal is the same, just presented at three different
%   levels. The purpose is to visualize the broading of the auditory
%   filters at higher input levels.
%
%   FIGURE 1 Greasy at 50 dB
%
%     This figure shows the DRNL of an input speech signal with a level
%     of 50 dB SPL.
%
%   FIGURE 2 Greasy at 70 dB
%
%     This figure shows the DRNL of an input speech signal with a level
%     of 50 dB SPL.
%
%   FIGURE 3 Greasy at 90 dB
%
%     This figure shows the DRNL of an input speech signal with a level
%     of 50 dB SPL.
%
%   See also: drnl, jepsen2008preproc


fs=16000;
siglen=7000;
% Increase the length of the test signal to add room for the filter delay
insig=postpad(greasy,siglen);

% The three levels in dB
lvl1=50;
lvl2=70;
lvl3=90;

% Set the dynamic range for plotting.
dynrange=80;

%% DRNL 
% Change 'bwmul' to generate more channels, this make the visualization
% pretty.
[outsig1, fc1] = drnl(setdbspl(insig,lvl1), fs, 'jepsen2008','bwmul',.1);
[outsig2, fc2] = drnl(setdbspl(insig,lvl2), fs, 'jepsen2008','bwmul',.1);
[outsig3, fc3] = drnl(setdbspl(insig,lvl3), fs, 'jepsen2008','bwmul',.1);

%% 'haircell' envelope extraction
outsig1 = ihcenvelope(outsig1,fs,'dau');
outsig2 = ihcenvelope(outsig2,fs,'dau');
outsig3 = ihcenvelope(outsig3,fs,'dau');

%% Lowpass filter the output for visualization and convert to dB
[mlp_b,mlp_a] = butter(2,50/(fs/2));
outsig1 = abs(filter(mlp_b,mlp_a,outsig1));
outsig2 = abs(filter(mlp_b,mlp_a,outsig2));
outsig3 = abs(filter(mlp_b,mlp_a,outsig3));

%% Visualization
ytick=[0,100,250,500,1000,2000,4000,8000];
ytickpos=freqtoerb(ytick);
xr=(0:siglen-1)/fs;
yr=linspace(freqtoerb(fc1(1)),freqtoerb(fc1(end)),length(fc1));
max1=20*log10(max(outsig1(:)));
max2=20*log10(max(outsig2(:)));
max3=20*log10(max(outsig3(:)));


figure(1);
imagesc(xr,yr,20*log10(outsig1).',[max1-dynrange,max1]); colorbar; axis('xy');
set(gca,'YTick',ytickpos);
set(gca,'YTickLabel',num2str(ytick(:)));

figure(2);
imagesc(xr,yr,20*log10(outsig2).',[max2-dynrange,max2]); colorbar; axis('xy');
set(gca,'YTick',ytickpos);
set(gca,'YTickLabel',num2str(ytick(:)));

figure(3);
imagesc(xr,yr,20*log10(outsig3).',[max3-dynrange,max3]); colorbar; axis('xy');
set(gca,'YTick',ytickpos);
set(gca,'YTickLabel',num2str(ytick(:)));

