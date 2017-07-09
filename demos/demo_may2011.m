%DEMO_MAY2011 Demo of the model estimating the azimuths of concurrent speakers
%
%
%   `demo_may2011` generates figures showing the result of the model estimating
%   the azimuth position of three concurrent speakers. Also, it returns the
%   estimated azimuths.
% 
%   Set `demo` to the following flags to shows other conditions:
%   
%   `1R`
%      one speaker in reverberant room
%
%   `2`
%      two speakers in free field
%
%   `3`
%      three speakers in free field (default)
%
%   `5`
%      five speakers in free field
%
%   .. figure::
%
%      Time-frequency-based azimuth estimates
%
%      This figure shows the azimuth estimates in the time-frequency
%      domain for three speakers.
%
%   .. figure::
%
%      Interaural time differences (ITDs)
%
%      This figure shows the ITDs in the time-frequency domain estimated
%      from the mixed signal of three concurrent speakers.
%
%   .. figure::
%
%      Interaural level differences (ILDs)
%
%      This figure shows the ILDs in the time-frequency domain estimated
%      from the mixed signal of three concurrent speakers.
%
%   .. figure::
%
%      Interaural coherence
%
%      This figure shows the interaural coherence in the time-frequency domain estimated
%      from the mixed signal of three concurrent speakers.
%
%   .. figure::
%
%      Frame-based azimuth estimates
%
%      This figure shows the azimuth directions in the time domain estimated
%      from the mixed signal of three concurrent speakers.
%
%   .. figure::
%
%      GMM pattern
%
%      This figure shows the pattern and the histogram obtained from the
%      GMM-estimator for the mixed signal of three concurrent speakers.
%
%   See also: may2011
%
%   References: may2011


%% Select binaural recordings
% 
% 
% Select one demo
if ~exist('demo','var')
  demo='3';
end
 
% Create signals
switch lower(demo)
    case '1r'
        % Create input signal
        [signal,fs] = sig_competingtalkers('one_speaker_reverb');
        
        % Find all active sources
        nSources = 1;
        
    case '2'
        % Create input signal
        [signal,fs] = sig_competingtalkers('two_of_three');
        
        % Find all active sources
        nSources = 2;
        
    case '3'
        % Create input signal
        [signal,fs] = sig_competingtalkers('three_of_three');
        
        % Find all active sources
        nSources = 3;
        
    case '5'
        % Create input signal
        [signal,fs] = sig_competingtalkers('five_speakers');
        
        % Find all active sources
        nSources = 5;
end


%% Perform GMM-based sound source localization
% 
% 
% Perform localization
out = may2011(signal,fs);


%% Plot results
% 
% 
% Plot time-frequency-based azimuth estimates
figure;
imagesc(out.azimuth,[-90 90]);
xlabel('Number of frames')
ylabel('Number of gammatone channels')
title('Time-frequency-based azimuth estimates')
colorbar;
may2011cbarlabel('Azimuth (deg)')
axis xy;

% Plot binaural cues
figure;
imagesc(out.itd,[-1e-3 1e-3]);
xlabel('Number of frames')
ylabel('Number of gammatone channels')
title('Interaural time difference (ITD)')
colorbar;
may2011cbarlabel('ITD (ms)')
axis xy;

figure;
imagesc(out.ild);
xlabel('Number of frames')
ylabel('Number of gammatone channels')
title('Interaural level difference (ILD)')
colorbar;
may2011cbarlabel('ILD (dB)')
axis xy;

figure;
imagesc(out.ic);
xlabel('Number of frames')
ylabel('Number of gammatone channels')
title('Interaural coherence (IC)')
colorbar;
may2011cbarlabel('IC')
axis xy;

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

% Histogram analysis of frame-based localization estimates
azEst=may2011estAzimuth_GMM(out,'HIST',nSources,1)