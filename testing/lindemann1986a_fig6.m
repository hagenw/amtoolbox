function lindemann1986a_fig6()
%LINDEMANN1986a_FIG6 Reproduces fig. 6 from lindemann1986a
%
%   LINDEMANN1986a_FIG6() reproduces fig.6 from lindemann1986a. Therefore the
%   cross-correlation of pure tone sinusoids with different ITDs is calculated. 
%   Because of the stationary character of the input signals T_int = inf is used
%   to produce only one time step in the crosscorr output from lindemann.
%   This is calculated for different ITDs and different inhibition factors c_s
%   (0,0.3,1). Afterwards for every c_s the correlation is plotted for every
%   used ITD dependend on the correlation-time delay.
%
%R lindemann1986a
%

%   AUTHOR: Hagen Wierstorf


% ------- Computation ----------------------------------------------------
% Sampling rate
fs = 44100;
% Frequency of the sinusoid
f = 500;

% Model parameter
T_int = inf;
w_f = 0;
M_f = 6; % not used, if w_f==0

% Calculate crosscorrelations for 21 ITD points between 0~ms and 1~ms
ii = 1;
cc_a = zeros(21,45);
cc_b = cc_a;
cc_c = cc_a;
for itd = linspace(0,1,21); 
    % Generate ITD shifted sinusoid
    sig = itdsin(f,itd,fs);
    % Use only the beginning of the signal to generate only one time instance of
    % the cross-correlation
    % FIXME: also we use only stationary inhibition here by varying c_s the
    % results depends much on the length of the used signal. I think Lindemann 
    % has not stated this fact nor calculated the fading time, so we have to do
    % this.
    sig = sig(1:ceil(0.01*fs),:);
    % Calculate cross-correlation for different inhibition factor c_s
    % (mean above all frequency channels)
    % fig. 6 (a)
    c_s = 0;
    cc_a(ii,:) = mean(lindemann(sig,fs,c_s,w_f,M_f,T_int),3);
    % fig. 6 (b)
    c_s = 0.3;
    cc_b(ii,:) = mean(lindemann(sig,fs,c_s,w_f,M_f,T_int),3);
    % fig. 6 (c)
    c_s = 1;
    cc_c(ii,:) = mean(lindemann(sig,fs,c_s,w_f,M_f,T_int),3);
    
    ii = ii+1;
end


% ------ Plotting --------------------------------------------------------
% Generate time axis
t = linspace(1,0,21);
tau = linspace(-1,1,45);
%
figure;
mesh(tau,t,cc_a);
view(0,57);
xlabel('correlation-time tau (ms)');
ylabel('interaural time difference (ms)');
tics('y',[0,0.2,0.4,0.6,0.8,1],['1';'0.8';'0.6';'0.4';'0.2';'0']);
tstr = sprintf('c_s = 0\nw_f = 0\nf = 500 Hz\n');
title(tstr);
%
figure;
mesh(tau,t,cc_b);
view(0,57);
xlabel('correlation-time tau (ms)');
ylabel('interaural time difference (ms)');
tics('y',[0,0.2,0.4,0.6,0.8,1],['1';'0.8';'0.6';'0.4';'0.2';'0']);
tstr = sprintf('c_s = 0.3\nw_f = 0\nf = 500 Hz\n');
title(tstr);
%
figure;
mesh(tau,t,cc_c);
view(0,57);
xlabel('correlation-time tau (ms)');
ylabel('interaural time difference (ms)');
tics('y',[0,0.2,0.4,0.6,0.8,1],['1';'0.8';'0.6';'0.4';'0.2';'0']);
tstr = sprintf('c_s = 1\nw_f = 0\nf = 500 Hz\n');
title(tstr);

