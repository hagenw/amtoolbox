function cc = lindemann1986a_fig11()
%LINDEMANN1986a_FIG Reproduces fig. 11 from lindemann1986a
%   Usage: cc = lindemann1986a_fig11;
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
fc = round(freqtoerb(f));   % corresponding frequency channel

% Model parameter
T_int = inf;
w_f = [0,0,0.035];
M_f = 6; % not used, if w_f==0
c_s = [0.3,1,0.3];

% Calculate crosscorrelations for 26 ILD points between 0~dB and 25~dB
nilds = 26; % number of used ILDs
ndl = 45;   % length of the delay line
ild = linspace(0,25,nilds);
cen = zeros(length(c_s),nilds);
for ii = 1:nilds 
    % Generate sinusoid with given ILD
    sig = ildsin(f,ild(ii),fs);
    % Use only the beginning of the signal to generate only one time instance of
    % the cross-correlation
    % NOTE: the signal has to be longer than N_1, which is 0.4*fs in this case
    % (see lindemann1986a p. 1614)
    sig = sig(1:ceil(0.41*fs),:);
    % Generate signal for normalization calculation
    if max(sig(:,1))>max(sig(:,2))
        sig0 = [sig(:,1) sig(:,1)];
    else
        sig0 = [sig(:,2) sig(:,2)];
    end
    % Calculate Psi_0 for normalization (see lindemann1986a p. 1614, eq. 25)
    tmp = squeeze(lindemann(sig,fs,0,0,M_f,T_int));
    % Store the needed frequency channel. NOTE: the cross-correlation
    % calculation starts with channel 5, so we have to subtract 5.
    cc_0 = tmp(:,fc-5);
    cc_0(ceil(ndl/2))
    % Calculate cross-correlation for different inhibition factor c_s 
    for jj = 1:length(c_s)
        % Calculate cross-correlation (and squeeze due to T_int==inf)
        tmp = squeeze(lindemann(sig,fs,c_s(jj),w_f(jj),M_f,T_int));
        % Store the needed frequency channel and normalize it
        cc = tmp(:,fc-5)/cc_0(ceil(ndl/2));
        % Calculate the position of the centroid
        cen(jj,ii) = centroid(cc);
    end
end


% ------ Plotting --------------------------------------------------------
% Generate time axis
tau = linspace(-1,1,ndl);
figure;
% Plot line for every condition
for jj = 1:length(c_s) 
    plot(ild,cen(jj,:));
    hold on;
end
xlabel('interaural level difference (dB)');
ylabel('displacement of the centroid d');
%tstr = sprintf('c_{inh} = %.1f\nw_f = 0.035\nf = 500 Hz\n',c_s(jj));
%title(tstr);

