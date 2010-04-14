function cen = lindemann1986a_fig7()
%LINDEMANN1986a_FIG7 Reproduces fig. 7 from lindemann1986a
%   Usage: cen = lindemann1986a_fig7;
%
%   Output parameters:
%       cen - Displacement of the centroid for different c_s values.
%             Dim: c_s x nitds
%
%   LINDEMANN1986a_FIG7() reproduces fig.7 from lindemann1986a. Therefore the
%   cross-correlation of pure tone sinusoids with different ITDs is calculated. 
%   Because of the stationary character of the input signals T_int = inf is used
%   to produce only one time step in the crosscorr output from lindemann.
%   This is calculated for different ITDs and different inhibition factors c_s
%   (0:0.2:1). Afterwards for every c_s the centroid of the auditory image is
%   calculated and plotted dependend on the ITD.
%
%   See also: lindemann, itdsin
%
%R lindemann1986a
%

%   AUTHOR: Hagen Wierstorf


% ------- Computation ----------------------------------------------------
% Sampling rate
fs = 44100;
% Frequency of the sinusoid
f = 500;
T = 1/f;
fc = round(freqtoerb(f));   % corresponding frequency channel

% Model parameter
T_int = inf;
w_f = 0;
M_f = 6; % not used, if w_f==0
c_s = 0:0.2:1;

% NOTE: the longer the signal, the more time we need for computation. On the
% other side N_1 needs to be long enough to eliminate any onset effects.
% Lindemann uses N_1 = 17640. Here I uses only N_1 = 2205 which gives the same
% results for this demo.
N_1 = ceil(25*T*fs);
siglen = ceil(30*T*fs);

% Calculate crosscorrelations for 21 ITD points between 0~ms and 1~ms
nitds = 21; % number of used ITDs
ndl = 2*round(fs/2000)+1;   % length of the delay line (see bincorr.m)
itd = linspace(0,1,nitds);
for ii = 1:nitds 
    % Generate ITD shifted sinusoid
    sig = itdsin(f,itd(ii),fs);
    % Use only the beginning of the signal to generate only one time instance of
    % the cross-correlation and apply an onset window
    sig = sig(1:siglen,:);
    sig = lindemannwin(sig,N_1);
    % Calculate cross-correlation for different inhibition factor c_s 
    for jj = 1:length(c_s)
        % Calculate cross-correlation (and squeeze due to T_int==inf)
        tmp = squeeze(lindemann(sig,fs,c_s(jj),w_f,M_f,T_int,N_1));
        % Store the needed frequency channel. NOTE: the cross-correlation
        % calculation starts with channel 5, so we have to subtract 4.
        cc = tmp(:,fc-4);
        % Calculate the position of the centroid
        cen(jj,ii) = centroid(cc);
    end
end


% ------ Plotting --------------------------------------------------------
figure;
for jj = 1:length(c_s)
    plot(itd,cen(jj,:));
    hold on;
end
xlabel('interaural time difference (ms)');
ylabel('displacement of the centroid d');
tstr = sprintf('w_f = 0\nf = %i Hz\n',f);
title(tstr);

