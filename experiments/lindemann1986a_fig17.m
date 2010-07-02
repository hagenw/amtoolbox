function [rcen rmax] = lindemann1986a_fig17()
%LINDEMANN1986a_FIG17 Reproduces fig. 17 from lindemann1986a
%   Usage: [rcen rmax] = lindemann1986a_fig17()
%
%   Output parameters:
%       rcen    - ITD value for getting the same lateralization with an ITD only
%                 stimulus compared to a stimulus with both ITD and ILD using
%                 the centroid of the cross-correlation.
%                 Dim: number of ITDs x number of ILDs
%       rmax    - ITD value for getting the same lateralization with an ITD only
%                 stimulus compared to a stimulus with both ITD and ILD using
%                 the maximum of the cross-correlation.
%                 Dim: number of ITDs x number of ILDs
%
%   LINDEMANN1986a_FIG17() reproduces fig.17 from lindemann1986a. Therefore the
%   cross-correlation of pure tone sinusoids with f=500 Hz with different ITDs
%   and with different combinations from ITDs and ILDs is calculated. 
%   Because of the stationary character of the input signals T_int = inf is used
%   to produce only one time step in the crosscorr output from lindemann.
%   After the calculation the values for the centroids and maxima of the ITD only
%   stimuli are searched to find the nearest value to the centroid and maxima of 
%   a given combined stimulus. The resulting ITD value is stored for different 
%   combinaition values.
%
%   See also: lindemann, itdsin, itdildsin
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
nitds_p = 21; % number of used ITDs for the ITD only stimuli
nitds_t = 4; % number of used ITDs for the combined stimuli
nilds_t = 21;  % number of used ILDs for the combined stimuli
ndl = 2*round(fs/2000)+1;     % length of the delay line (see bincorr.m)
itd_p = linspace(-0.36,0.72,nitds_p);
ild_t = linspace(-9,9,nilds_t);
itd_t = [0,0.09,0.18,0.27];

% Calculate the centroids for the ITD only stimuli
cen_p = zeros(1,nitds_p);
max_p = zeros(1,nitds_p);
for ii = 1:nitds_p
    % Generate sinusoid with given ILD
    sig = itdsin(f,itd_p(ii),fs);
    % Use only the beginning of the signal to generate only one time instance of
    % the cross-correlation and apply a linear onset window
    sig = sig(1:siglen,:);
    sig = lindemannwin(sig,N_1);
    % Calculate cross-correlation (and squeeze due to T_int==inf)
    tmp = squeeze(lindemann(sig,fs,c_s,w_f,M_f,T_int,N_1));
    % Store the needed frequency channel
    cc = tmp(:,fc-4);
    % Find the maximum position
    max_p(ii) = findmax(cc);
    % Calculate the position of the centroid
    cen_p(ii) = lindcentroid(cc);
end

% Calculate the centroids for the combined stimuli
cen_t = zeros(nitds_t,nilds_t);
max_t = zeros(nitds_t,nilds_t);
for ii = 1:nilds_t
    for jj = 1:nitds_t
        sig = itdildsin(f,itd_t(jj),ild_t(ii),fs);
        sig = sig(1:siglen,:);
        sig = lindemannwin(sig,N_1);
        tmp = squeeze(lindemann(sig,fs,c_s,w_f,M_f,T_int,N_1));
        cc = tmp(:,fc-4);
        max_t(jj,ii) = findmax(cc);
        cen_t(jj,ii) = lindcentroid(cc);
    end
end

% ------ Fiting ---------------------------------------------------------
% For the results for the combined centroids find the nearest centroid for the
% ITD only stimuli
r = zeros(nitds_t,nilds_t);
for ii = 1:nilds_t
    for jj = 1:nitds_t
        idx = findnearest(cen_p,cen_t(jj,ii));
        rcen(jj,ii) = itd_p(idx);
        idx = findnearest(max_p,max_t(jj,ii));
        rmax(jj,ii) = itd_p(idx);
    end
end


% ------ Plotting --------------------------------------------------------
% First plot the only data from experiments
figure; % fig 17 (a)
d = lindemann1986a_data(17);
plot(d(:,1),d(:,5),'x-r', ...   %  0 ms
     d(:,1),d(:,4),'x-b', ...   %  0.09 ms
     d(:,1),d(:,3),'x-g', ...   %  0.18 ms
     d(:,1),d(:,2),'x-b')       %  0.27 ms
legend('0.27 ms','0.18 ms','0.09 ms','0 ms');
axis([-9 9 -0.36 0.72]);
set(gca,'XTick',-9:3:9);
xlabel('interaural level difference (dB)');
ylabel('interaural time difference (ms)');
% Then plot model results
figure; % fig 17 (b)
% Plot line for every condition
for jj = 1:nitds_t
    plot(ild_t,rcen(jj,:));
    hold on;
end
axis([-9 9 -0.36 0.72]);
set(gca,'XTick',-9:3:9);
xlabel('interaural level difference (dB)');
ylabel('interaural time difference (ms)');
%
figure; % fig 17 (c)
for jj = 1:nitds_t
    plot(ild_t,rmax(jj,:));
    hold on;
end
axis([-9 9 -0.36 0.72]);
set(gca,'XTick',-9:3:9);
xlabel('interaural level difference (dB)');
ylabel('interaural time difference (ms)');

% ------ Subfunctions ----------------------------------------------------
function idx = findnearest(list,val)
[dif,idx] = min(abs(list-val));

function m = findmax(list)
[m,idx] = max(list);
d = linspace(-1,1,length(list));
m = d(idx);
