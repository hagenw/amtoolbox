%% Parameter Settings

% model fiber parameters
CF    = 1.5e3;   % CF in Hz;   
fiberType = 1;  % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High

% stimulus parameters
F0 = CF;     % stimulus frequency in Hz
fsstim = 100e3;  % sampling rate in Hz (must be 100, 200 or 500 kHz)
T  = 50e-3;  % stimulus duration in seconds
rt = 2.5e-3; % rise/fall time in seconds
stimdb = 65; % stimulus intensity in dB SPL

% peri-stimulus time histogram (PSTH) parameters
nrep = 1;               % number of stimulus repetitions (e.g., 50);
psthbinwidth = 0.5e-3;  % binwidth in seconds;


%% Computations

% Stimulus generation
t = 0:1/fsstim:T-1/fsstim; % time vector
mxpts = length(t);
irpts = rt*fsstim;
stim = sqrt(2)*20e-6*10^(stimdb/20)*sin(2*pi*F0*t); % unramped stimulus
stim(1:irpts)= stim(1:irpts).*(0:(irpts-1))/irpts; 
stim((mxpts-irpts):mxpts)=stim((mxpts-irpts):mxpts).*(irpts:-1:0)/irpts;

% AN modeling
[ANresp,fc,vihc,psth] = zilany2013(...
  stimdb,stim,fsstim,...
  'flow',CF','fhigh',CF,'nfibers',1,'fiberType',fiberType);

% PSTH conversion
timeout = (1:length(psth))*1/fsstim;
psthbins = round(psthbinwidth*fsstim);  % number of psth bins per psth bin
psthtime = timeout(1:psthbins:end); % time vector for psth
pr = sum(reshape(psth,psthbins,length(psth)/psthbins))/nrep; % pr of spike in each bin
Psth = pr/psthbinwidth; % psth in units of spikes/s


%% Plots

figure
subplot(4,1,1)
plot(timeout,[stim zeros(1,length(timeout)-length(stim))])
title('Input Stimulus')
ylabel('Pascal')

subplot(4,1,2)
plot(timeout,vihc(1:length(timeout)))
title('IHC Output')
ylabel('Volts')

subplot(4,1,3)
plot(timeout,ANresp);
xl = xlim;
title('Mean Rate Output')
ylabel('spikes/s')

subplot(4,1,4)
bar(psthtime,Psth)
xlim(xl)
title('Peri-stimulus Time Histogram')
xlabel('Time (s)')
ylabel('spikes/s')
