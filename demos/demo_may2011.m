%DEMO_MAY2011 Demo of the May et al. (2011)
%  
%   This script does not work unless speakers1234 is provided
%
%   %   .. figure::
%
%   See also: may2011


% Load example

% load speakers24
% load speakers123
load speakers1234
% load speakers12345

% Sampling frequency
fs = 16000;

% Perform localization
out = may2011(signal,fs);


%% Plot results
% 
% 
% Plot frame-based azimuth estimates
figure;
plot(out.azFrames,'k.','linewidth',2);
xlabel('Number of frames')
ylabel('Azimuth (deg)')
title('Frame-based azimuth estimates')
xlim([-inf inf])
ylim([-90 90])
grid on;
axis xy;

% Find all active sources, make no assumptions
nSources = inf;

% Histogram analysis of frame-based localization estimates
estAzimuth_GMM(out,'HIST',nSources);

% Plot time-frequency-based azimuth estimates
figure;
imagesc(out.azimuth,[-90 90]);
xlabel('Number of frames')
ylabel('Number of gammatone channels')
title('Time-frequency-based azimuth estimates')
colorbar;
cbarlabel('Azimuth (deg)')
axis xy;

% Plot binaural cues
figure;
imagesc(out.itd,[-1e-3 1e-3]);
xlabel('Number of frames')
ylabel('Number of gammatone channels')
title('Interaural time difference (ITD)')
colorbar;
cbarlabel('ITD (ms)')
axis xy;

figure;
imagesc(out.ild);
xlabel('Number of frames')
ylabel('Number of gammatone channels')
title('Interaural level difference (ILD)')
colorbar;
cbarlabel('ILD (dB)')
axis xy;

figure;
imagesc(out.ic);
xlabel('Number of frames')
ylabel('Number of gammatone channels')
title('Interaural coherence (IC)')
colorbar;
cbarlabel('IC')
axis xy;

% Arrange open figures on screen
arrange(1:6,[3 2])