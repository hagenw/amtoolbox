function cen = lindemann1986a_fig14b()
%LINDEMANN1986a_FIG14b Reproduces fig. 14 (b) from lindemann1986a
%   Usage: cen = lindemann1986a_fig14b()
%
%   Output parameters:
%       cen - Displacement of the centroid as a function of the ITD averaged
%             over different small ILD with a standard deviation of 0,1,2,3,4,5.
%             Dim: nilds x nitds
%
%   LINDEMANN1986a_FIG14b() reproduces fig.14 (b) from lindemann1986a. Therefore
%   the cross-correlation of pure tone sinusoids with f=500 Hz with different
%   ILDs and an ITD of 1~ms is calculated. Because of the stationary character
%   of the input signals T_int = inf is used to produce only one time step in
%   the crosscorr output from lindemann. This is calculated for different samll
%   ILDs with a standard deviation of 0-5. Afterwards for every standard 
%   deviation the mean centroid displacement is calculated and plotted
%   dependend on the ITD.
%
%   See also: lindemann, itdildsin
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
c_s = 1.0;

% NOTE: the longer the signal, the more time we need for computation. On the
% other side N_1 needs to be long enough to eliminate any onset effects.
% Lindemann uses N_1 = 17640. Here I uses only N_1 = 2205 which gives the same
% results for this demo.
N_1 = ceil(25*T*fs);
siglen = ceil(30*T*fs);

fprintf(1,'NOTE: this test function will need a lot of time!\n\n');

% Calculate crosscorrelations for 21 ITD points between 0~ms and 1~ms
nitds = 21; % number of used ITDS
nilds = 201; % number of used ILDs
ndl = 2*round(fs/2000)+1;   % length of the delay line (see bincorr.m)
itd = linspace(0,1,nitds);
ild_std = [1,2,3,4,5];
cen = zeros(length(ild_std)+1,nitds);
centmp = zeros(nilds,nitds);
for ii = 1:nitds
    % Show progress
    progressbar(ii,nitds);
    % First generate the result for std(ILD) == 0
    sig = itdsin(f,itd(ii),fs);
    sig = sig(1:siglen,:);
    sig = lindemannwin(sig,N_1);
    tmp = squeeze(lindemann(sig,fs,c_s,w_f,M_f,T_int,N_1));
    cc = tmp(:,fc-4);
    cen(1,ii) = centroid(cc);
    % Generate results for std(ILD) ~= 0
    for nn = 1:length(ild_std)
        % Generate normal distributed ILDs with mean ~ 0 and std = 1
        tmp = randn(nilds,1);
        tmp = tmp/std(tmp);
        % Generate desired ILD distribution
        ild = tmp * ild_std(nn);
        % For all distributed ILD values calculate the centroid
        for jj = 1:nilds
            sig = itdildsin(f,itd(ii),ild(jj),fs);
            sig = sig(1:siglen,:);
            sig = lindemannwin(sig,N_1);
            tmp = squeeze(lindemann(sig,fs,c_s,w_f,M_f,T_int,N_1));
            cc = tmp(:,fc-4);
            centmp(jj,ii) = centroid(cc);
        end
        % Calculate the mean centroid above the ILD distribution
        cen(nn+1,ii) = mean(centmp(:,ii));
    end
end


% ------ Ploting --------------------------------------------------------
figure
for nn = 1:length(ild_std)+1
    plot(itd,cen(nn,:)); hold on;
end
axis([0 1 0 1]);
xlabel('interaural level difference (dB)');
ylabel('displacement of the centroid d');
tstr = sprintf('w_f = 0\nc_{inh} = 1\nf = %i Hz\n',f);
title(tstr);

