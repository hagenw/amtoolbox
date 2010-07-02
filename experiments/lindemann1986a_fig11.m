function cc = lindemann1986a_fig11()
%LINDEMANN1986a_FIG11 Reproduces fig. 11 from lindemann1986a
%   Usage: cc = lindemann1986a_fig11()
%
%   Output parameters:
%       cc  - cross-correlation result of the to figure.
%             Dim: number of c_s conditions x nilds x delay line length
%
%   LINDEMANN1986a_FIG11() reproduces fig.11 from lindemann1986a. Therefore the
%   cross-correlation of pure tone sinusoids with different ILDs is calculated. 
%   Because of the stationary character of the input signals T_int = inf is used
%   to produce only one time step in the crosscorr output from lindemann.
%   This is calculated for different ILDs, an inhibition factor c_s = 0.3 and a
%   monaural detector factor w_f = 0.035. Afterwards for every c_s the ILD is
%   plotted dependend on the correaltion time.
%
%   See also: lindemann, ildsin
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
w_f = [0,0,0.035];
M_f = 6; % not used, if w_f==0
c_s = [0.3,1,0.3];

% NOTE: the longer the signal, the more time we need for computation. On the
% other side N_1 needs to be long enough to eliminate any onset effects.
% Lindemann uses N_1 = 17640. Here I uses only N_1 = 2205 which gives the same
% results for this demo.
N_1 = ceil(25*T*fs);
siglen = ceil(30*T*fs);

% Calculate crosscorrelations for 26 ILD points between 0~dB and 25~dB
nilds = 26; % number of used ILDs
ndl = 2*round(fs/2000)+1;   % length of the delay line (see bincorr.m)
ild = linspace(0,25,nilds);
cen = zeros(length(c_s),nilds);
for ii = 1:nilds 
    % Generate sinusoid with given ILD
    sig = ildsin(f,ild(ii),fs);
    % Use only the beginning of the signal to generate only one time instance of
    % the cross-correlation and apply onset window
    sig = sig(1:siglen,:);
    sig = lindemannwin(sig,N_1);
    % Calculate cross-correlation for different inhibition factor c_s 
    for jj = 1:length(c_s)
        % Calculate cross-correlation (and squeeze due to T_int==inf)
        tmp = squeeze(lindemann(sig,fs,c_s(jj),w_f(jj),M_f,T_int,N_1));
        % Store the needed frequency channel
        cc = tmp(:,fc-4);
        % Calculate the position of the centroid
        cen(jj,ii) = lindcentroid(cc);
    end
end


% ------ Plotting --------------------------------------------------------
% Generate time axis
tau = linspace(-1,1,ndl);
figure;
% Plot data from experiments
data = lindemann1986a_data(11,'yost');
plot(data(:,1),data(:,2)*0.058,'+'); hold on;
data = lindemann1986a_data(11,'sayers');
plot(data(:,1),data(:,2)*0.058,'*');
% Plot line for every condition
for jj = 1:length(c_s) 
    plot(ild,cen(jj,:));
end
axis([0 25 0 0.8]);
legend('500 Hz, Yost','600 Hz, Sayers','500 Hz, model');
xlabel('interaural level difference (dB)');
ylabel('displacement of the centroid d');

