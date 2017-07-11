%demo_lopezpoveda2001  Widening of filters in the lopezpoveda2001
%
%   This script displays three auditory spectrograms generated using the
%   lopezpoveda2001. The input signal is the same, just presented at three different
%   levels. The purpose is to visualize the broading of the auditory
%   filters at higher input levels.
%
%   .. figure::
%
%      Greasy at 50 dB
%
%      This figure shows the lopezpoveda2001 of an input speech signal with a level
%      of 50 dB SPL.
%
%   .. figure::
%
%      Greasy at 70 dB
%
%      This figure shows the lopezpoveda2001 of an input speech signal with a level
%      of 50 dB SPL.
%
%   .. figure::
%
%      Greasy at 90 dB
%
%      This figure shows the lopezpoveda2001 of an input speech signal with a level
%      of 50 dB SPL.
%
%   See also: lopezpoveda2001


fs=32000;
siglen=7000;
% Increase the length of the test signal to add room for the filter delay
insig=postpad(greasy,siglen);

% The three levels in dB
lvl1=50;
lvl2=70;
lvl3=90;

% Set the dynamic range for plotting.
dynrange=30;

%% lopezpoveda2001 
% Change 'bwmul' to generate more channels, this make the visualization
% pretty.
[outsig1, fc1] = lopezpoveda2001(setdbspl(insig,lvl1), fs, 'bwmul',.1);
[outsig2, fc2] = lopezpoveda2001(setdbspl(insig,lvl2), fs, 'bwmul',.1);
[outsig3, fc3] = lopezpoveda2001(setdbspl(insig,lvl3), fs, 'bwmul',.1);

%% 'haircell' envelope extraction
outsig1 = ihcenvelope(outsig1,fs,'ihc_dau');
outsig2 = ihcenvelope(outsig2,fs,'ihc_dau');
outsig3 = ihcenvelope(outsig3,fs,'ihc_dau');

%% Lowpass filter the output for visualization and convert to dB
[mlp_b,mlp_a] = butter(2,50/(fs/2));
outsig1 = abs(filter(mlp_b,mlp_a,outsig1));
outsig2 = abs(filter(mlp_b,mlp_a,outsig2));
outsig3 = abs(filter(mlp_b,mlp_a,outsig3));

%% Visualization
ftz=20;

figure(1);
plotfilterbank(outsig1,1,fc1,fs,dynrange,'audtick');
xlabel('Time (s)','fontsize',ftz,'fontname','arial')
ylabel('Frequency (Hz)','fontsize',ftz)
colormap(flipud(bone));

figure(2);
plotfilterbank(outsig2,1,fc1,fs,dynrange,'audtick');
colormap(flipud(bone));

figure(3);
plotfilterbank(outsig3,1,fc1,fs,dynrange,'audtick');
colormap(flipud(bone));