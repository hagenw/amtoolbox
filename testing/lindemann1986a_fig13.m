function r = lindemann1986a_fig13()
%LINDEMANN1986a_FIG13 Reproduces fig. 13 from lindemann1986a
%   Usage: r = lindemann1986a_fig13()
%
%   Output parameters:
%       r   - ILD value for getting the same lateralization with an ILD only
%             stimulus compared to a stimulus with both ITD and ILD.
%             Dim: number of ILDs x number of ITDs
%
%   LINDEMANN1986a_FIG13() reproduces fig.13 from lindemann1986a. Therefore the
%   cross-correlation of pure tone sinusoids with f=500 Hz with different ILDs
%   and with different combinations from ITDs and ILDs is calculated. 
%   Because of the stationary character of the input signals T_int = inf is used
%   to produce only one time step in the crosscorr output from lindemann.
%   After the calculation the values for the centroids of the ILD only stimuli
%   are searched to find the nearest value to the centroid of a given combined
%   stimulus. The resulting ILD value is stored for different combinaition
%   values.
%
%   See also: lindemann, ildsin, itdildsin
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
nilds_p = 21; % number of used ILDs for the ILD only stimuli
nitds_t = 21; % number of used ITDs for the combined stimuli
nilds_t = 6;  % number of used ILDs for the combined stimuli
ndl = 2*round(fs/2000)+1;     % length of the delay line (see bincorr.m)
ild_p = linspace(-10,40,nilds_p);
itd_t = linspace(-1,1,nitds_t);
ild_t = [-3,0,3,9,15,25];

% Calculate the centroids for the ILD only stimuli
cen_p = zeros(1,nilds_p);
for ii = 1:nilds_p
    % Generate sinusoid with given ILD
    sig = ildsin(f,ild_p(ii),fs);
    % Use only the beginning of the signal to generate only one time instance of
    % the cross-correlation and apply a linear onset window
    sig = sig(1:siglen,:);
    sig = lindemannwin(sig,N_1);
    % Calculate cross-correlation (and squeeze due to T_int==inf)
    tmp = squeeze(lindemann(sig,fs,c_s,w_f,M_f,T_int,N_1));
    % Store the needed frequency channel
    cc = tmp(:,fc-4);
    % Calculate the position of the centroid
    cen_p(ii) = centroid(cc);
end

% Calculate the centroids for the combined stimuli
cen_t = zeros(nilds_t,nitds_t);
for ii = 1:nitds_t
    for jj = 1:nilds_t
        sig = itdildsin(f,itd_t(ii),ild_t(jj),fs);
        sig = sig(1:siglen,:);
        sig = lindemannwin(sig,N_1);
        tmp = squeeze(lindemann(sig,fs,c_s,w_f,M_f,T_int,N_1));
        cc = tmp(:,fc-4);
        cen_t(jj,ii) = centroid(cc);
    end
end

% ------ Fiting ---------------------------------------------------------
% For the results for the combined centroids find the nearest centroid for the
% ILD only stimuli
r = zeros(nilds_t,nitds_t);
for ii = 1:nitds_t
    for jj = 1:nilds_t
        idx = findnearest(cen_p,cen_t(jj,ii));
        r(jj,ii) = ild_p(idx);
    end
end


% ------ Plotting --------------------------------------------------------
% Generate time axis
tau = linspace(-1,1,ndl);
% First plot the only data from experiments
figure;
d = lindemann1986a_data(13);
plot(d(:,1),d(:,2),'x-r', ...   % -3dB
     d(:,1),d(:,3),'x-b', ...   %  3dB
     d(:,1),d(:,4),'x-g', ...   %  9dB
     d(:,1),d(:,5),'x-b', ...   % 15dB
     d(:,1),d(:,6),'x-r')       % 25dB
legend('25dB','15dB','9dB','3dB','-3dB');
axis([-1 1 -10 40]);
set(gca,'XTick',-1:0.4:1);
xlabel('interaural time difference (ms)');
ylabel('interaural level difference (dB)');
% Then plot model results
figure;
% Plot line for every condition
for jj = 1:nilds_t
    plot(itd_t,r(jj,:));
    hold on;
end
axis([-1 1 -10 40]);
set(gca,'XTick',-1:0.4:1);
xlabel('interaural time difference (ms)');
ylabel('interaural level difference (dB)');


% ------ Subfunctions ----------------------------------------------------
function idx = findnearest(list,val)
[dif,idx] = min(abs(list-val));
