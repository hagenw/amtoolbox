function cc = lindemann1986a_fig18()
%LINDEMANN1986a_FIG18 Reproduces fig. 18 from lindemann1986a
%   Usage: cc = lindemann1986a_fig18()
%
%   Output parameters:
%       cc  - cross-correlation result of the figure.
%             Dim: number of interaural coherences x delay line length
%
%   LINDEMANN1986a_FIG18() reproduces fig.18 from lindemann1986a. Therefore the
%   cross-correlation of pink noise with different interaural coherence values
%   is calculated. 
%   Because of the stationary character of the input signals T_int = inf is used
%   to produce only one time step in the crosscorr output from lindemann.
%   Afterwards for every interaural coherence value the correlation is plotted 
%   dependend on the correlation-time delay.
%
%   See also: lindemann, corpinknoise
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
w_f = 0.035;
M_f = 6; % not used, if w_f==0
c_s = 0.11;

% NOTE: the longer the signal, the more time we need for computation. On the
% other side N_1 needs to be long enough to eliminate any onset effects.
% Lindemann uses N_1 = 17640. Here I uses only N_1 = 2205 which gives the same
% results for this demo.
N_1 = ceil(25*T*fs);
siglen = ceil(30*T*fs);

% Calculate crosscorrelations for 21 ITD points between 0~ms and 1~ms
niacs = 11; % number of used interaural coherence values
ndl = 2*round(fs/2000)+1;   % length of the delay line (see bincorr.m)
iac = linspace(0,1,niacs);
cc = zeros(niacs,ndl);
for ii = 1:niacs; 
    % Generate ITD shifted sinusoid
    sig = corpinknoise(fs,iac(ii));
    % Aplly onset window
    sig = lindemannwin(sig,N_1);
    % Calculate cross-correlation (and squeeze due to T_int==inf)
    tmp = squeeze(lindemann(sig,fs,c_s,w_f,M_f,T_int,N_1));
    % Store the needed frequency channel. NOTE: the cross-correlation
    % calculation starts with channel 5, so we have to subtract 4.
    cc(ii,:) =  tmp(:,fc-4);
end


% ------ Plotting --------------------------------------------------------
% Generate time axis
tau = linspace(-1,1,ndl);
% Plot figure for every c_s condition
figure;
mesh(tau,iac,cc);
view(0,57);
xlabel('correlation-time tau (ms)');
ylabel('degree of interaural coherence');

