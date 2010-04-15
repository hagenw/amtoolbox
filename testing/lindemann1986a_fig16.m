function r = lindemann1986a_fig16()
%LINDEMANN1986a_FIG16 Reproduces fig. 16 from lindemann1986a
%   Usage: cc = lindemann1986a_fig16()
%
%   Output parameters:
%       r   - The ILD value for which the cross-correlation for a combined ITD,
%             ILD stimulus has two peaks with the same (nearest) height,
%             dependend on the ITD
%             Dim: number of ITDs x 1
%
%   LINDEMANN1986a_FIG16() reproduces fig.16 from lindemann1986a. Therefore the
%   cross-correlation of pure tone sinusoids with f=500 Hz with different 
%   combinations from ITDs and ILDs is calculated. 
%   Because of the stationary character of the input signals T_int = inf is used
%   to produce only one time step in the crosscorr output from lindemann.
%   After the calculation the combinations of ILD and ITD are looked for that
%   have two peaks in its cross-correlation which have the nearly same height.
%   The corresponding ILD value is than stored in r dependend on the ITD.
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
N_1 = 1;
w_f = 0.035;
M_f = 6; % not used, if w_f==0
c_s = 0.3;

% NOTE: the longer the signal, the more time we need for computation. On the
% other side N_1 needs to be long enough to eliminate any onset effects.
% Lindemann uses N_1 = 17640. Here I uses only N_1 = 2205 which gives the same
% results for this demo.
N_1 = ceil(25*T*fs);
siglen = ceil(30*T*fs);

% Calculate crosscorrelations for 26 ILD points between 0~dB and 25~dB
nilds = 21; % number of used ILDs for the ILD only stimuli
nitds = 21; % number of used ITDs for the combined stimuli
ndl = 2*round(fs/2000)+1;     % length of the delay line (see bincorr.m)
ild = linspace(0,10,nilds);
itd = linspace(-1,0,nitds);

% Calculate the centroids for ILD+ITD stimuli
cc = zeros(nitds,nilds,ndl);
for ii = 1:nitds
    for jj = 1:nilds
        % Generate sinusoid with given ILD
        sig = itdildsin(f,itd(ii),ild(jj),fs);
        % Use only the beginning of the signal to generate only one time 
        % instance of the cross-correlation and apply a linear onset window
        sig = sig(1:siglen,:);
        sig = lindemannwin(sig,N_1);
        % Calculate cross-correlation (and squeeze due to T_int==inf)
        tmp = squeeze(lindemann(sig,fs,c_s,w_f,M_f,T_int,N_1));
        % Store the needed frequency channel
        cc(ii,jj,:) = tmp(:,fc-4);
    end
end


% ------ Fiting ----------------------------------------------------------
% For every ITD find the ILD with gives a centroid near 0
r = zeros(nitds,1);
for ii = 1:nitds
    % Find two peaks of same height
    idx = findpeaks(squeeze(cc(ii,:,:)));
    r(ii) = ild(idx);
end


% ------ Plotting --------------------------------------------------------
figure;
data = lindemann1986a_data(16);
plot(data(:,1),data(:,2),'+'); hold on;
plot(itd,r);
axis([-1 0 0 15]);
xlabel('interaural time difference (ms)');
ylabel('interaural level difference (dB)');


% ------ Subfunctions ----------------------------------------------------
function idx = findpeaks(cc)
% FINDPEAKS Finds two peaks of the same height in a given cross-correlation
h = zeros(size(cc,1),1);
for ii = 1:size(cc,1)
    % Find a peak in the two halfes of the cross-correlation and subtract them
    h(ii) = max(cc(ii,1:30)) - max(cc(ii,31:end));
end
% Find the index with the smallest distance between the two peaks
[val,idx] = min(abs(h));
